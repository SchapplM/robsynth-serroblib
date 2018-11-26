% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S2RR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
%
% Output:
% JaD_rot [3x2]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S2RR1_jacobiaD_rot_2_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_jacobiaD_rot_2_floatb_twist_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_jacobiaD_rot_2_floatb_twist_sym_varpar: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_jacobiaD_rot_2_floatb_twist_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:44:42
% EndTime: 2018-11-16 16:44:42
% DurationCPUTime: 0.03s
% Computational Cost: add. (18->3), mult. (82->10), div. (42->6), fcn. (80->3), ass. (0->8)
t51 = sin(qJ(1));
t53 = cos(qJ(1));
t62 = 0.1e1 / t51 ^ 2 * t53 ^ 2;
t65 = -0.1e1 - t62;
t39 = cos(atan2(0, t51));
t38 = 0.1e1 / t39 ^ 2;
t33 = t38 * t62 + 0.1e1;
t1 = [0, 0; 0.2e1 * (0.1e1 / t33 + t65 / t33 ^ 2 * t38) * qJD(1) * t65 / t51 * t53 / t39, 0; 0, 0;];
JaD_rot  = t1;
