% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S3RPR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
%
% Output:
% JaD_rot [3x3]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S3RPR1_jacobiaD_rot_3_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_jacobiaD_rot_3_floatb_twist_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_jacobiaD_rot_3_floatb_twist_sym_varpar: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_jacobiaD_rot_3_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:34
% EndTime: 2018-11-14 10:14:34
% DurationCPUTime: 0.08s
% Computational Cost: add. (128->13), mult. (322->30), div. (40->4), fcn. (356->4), ass. (0->21)
t54 = sin(qJ(3));
t55 = cos(qJ(1));
t67 = sin(qJ(1));
t68 = cos(qJ(3));
t59 = t54 * t67 + t55 * t68;
t44 = 0.1e1 / t59 ^ 2;
t70 = t44 * t59;
t69 = qJD(1) - qJD(3);
t43 = 0.1e1 / t59;
t58 = -t54 * t55 + t67 * t68;
t38 = t69 * t58;
t42 = t58 ^ 2;
t64 = t42 * t44;
t41 = 0.1e1 + t64;
t65 = t69 * t70;
t62 = t58 * t65;
t45 = t43 * t44;
t63 = t42 * t45;
t66 = (t38 * t63 + t62) / t41 ^ 2;
t39 = 0.1e1 / t41;
t1 = [0, 0, 0; 0, 0, 0; 0.2e1 * (t43 * t59 + t64) * t66 + (-0.2e1 * t62 - (-t43 + 0.2e1 * t63 + t70) * t38) * t39, 0, -0.2e1 * t66 - 0.2e1 * (-t39 * t65 - (t38 * t39 * t45 - t44 * t66) * t58) * t58;];
JaD_rot  = t1;
