% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPPR2
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
%
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S4RPPR2_jacobiaD_rot_4_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_jacobiaD_rot_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_jacobiaD_rot_4_floatb_twist_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_jacobiaD_rot_4_floatb_twist_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:43
% EndTime: 2018-11-14 13:47:43
% DurationCPUTime: 0.09s
% Computational Cost: add. (306->14), mult. (322->30), div. (40->4), fcn. (356->4), ass. (0->22)
t70 = pkin(6) + qJ(4);
t60 = sin(t70);
t61 = cos(qJ(1));
t67 = cos(t70);
t75 = sin(qJ(1));
t65 = t75 * t60 + t61 * t67;
t50 = 0.1e1 / t65 ^ 2;
t77 = t50 * t65;
t76 = qJD(1) - qJD(4);
t49 = 0.1e1 / t65;
t64 = -t61 * t60 + t75 * t67;
t44 = t76 * t64;
t48 = t64 ^ 2;
t72 = t48 * t50;
t47 = 0.1e1 + t72;
t73 = t76 * t77;
t69 = t64 * t73;
t51 = t49 * t50;
t71 = t48 * t51;
t74 = (t44 * t71 + t69) / t47 ^ 2;
t45 = 0.1e1 / t47;
t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0.2e1 * (t49 * t65 + t72) * t74 + (-0.2e1 * t69 - (-t49 + 0.2e1 * t71 + t77) * t44) * t45, 0, 0, -0.2e1 * t74 - 0.2e1 * (-t45 * t73 - (t44 * t45 * t51 - t50 * t74) * t64) * t64;];
JaD_rot  = t1;
