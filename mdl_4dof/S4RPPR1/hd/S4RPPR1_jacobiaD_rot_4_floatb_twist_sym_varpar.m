% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPPR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
%
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S4RPPR1_jacobiaD_rot_4_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_jacobiaD_rot_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_jacobiaD_rot_4_floatb_twist_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_jacobiaD_rot_4_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:46:52
% EndTime: 2018-11-14 13:46:52
% DurationCPUTime: 0.09s
% Computational Cost: add. (306->14), mult. (322->30), div. (40->4), fcn. (356->4), ass. (0->22)
t67 = qJ(1) + pkin(6);
t57 = cos(t67);
t58 = sin(qJ(4));
t64 = sin(t67);
t72 = cos(qJ(4));
t62 = t57 * t72 + t64 * t58;
t47 = 0.1e1 / t62 ^ 2;
t74 = t47 * t62;
t73 = qJD(1) - qJD(4);
t46 = 0.1e1 / t62;
t61 = -t57 * t58 + t64 * t72;
t41 = t73 * t61;
t45 = t61 ^ 2;
t69 = t45 * t47;
t44 = 0.1e1 + t69;
t70 = t73 * t74;
t66 = t61 * t70;
t48 = t46 * t47;
t68 = t45 * t48;
t71 = (t41 * t68 + t66) / t44 ^ 2;
t42 = 0.1e1 / t44;
t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0.2e1 * (t46 * t62 + t69) * t71 + (-0.2e1 * t66 - (-t46 + 0.2e1 * t68 + t74) * t41) * t42, 0, 0, -0.2e1 * t71 - 0.2e1 * (-t42 * t70 - (t41 * t42 * t48 - t47 * t71) * t61) * t61;];
JaD_rot  = t1;
