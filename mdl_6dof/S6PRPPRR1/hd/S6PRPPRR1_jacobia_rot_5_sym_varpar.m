% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.08s
% Computational Cost: add. (252->19), mult. (626->46), div. (35->9), fcn. (880->13), ass. (0->35)
t65 = sin(pkin(10));
t66 = sin(pkin(6));
t75 = t65 * t66;
t69 = cos(pkin(6));
t64 = sin(pkin(11));
t67 = cos(pkin(11));
t70 = sin(qJ(2));
t71 = cos(qJ(2));
t73 = t71 * t64 + t70 * t67;
t58 = t73 * t69;
t59 = t70 * t64 - t71 * t67;
t68 = cos(pkin(10));
t53 = -t65 * t58 - t68 * t59;
t63 = pkin(12) + qJ(5);
t61 = sin(t63);
t62 = cos(t63);
t48 = t53 * t62 + t61 * t75;
t46 = 0.1e1 / t48 ^ 2;
t47 = t53 * t61 - t62 * t75;
t74 = t47 ^ 2 * t46 + 0.1e1;
t72 = t59 * t69;
t57 = t73 * t66;
t56 = t59 * t66;
t55 = 0.1e1 / t56 ^ 2;
t51 = t65 * t72 - t68 * t73;
t50 = -t65 * t73 - t68 * t72;
t49 = -t68 * t58 + t65 * t59;
t45 = atan2(t50, t56);
t43 = cos(t45);
t42 = sin(t45);
t41 = 0.1e1 / t74;
t40 = t42 * t50 + t43 * t56;
t39 = 0.1e1 / t40 ^ 2;
t37 = (t49 / t56 - t57 * t50 * t55) / (t50 ^ 2 * t55 + 0.1e1);
t1 = [0, t37, 0, 0, 0, 0; 0 (t53 / t40 + (t42 * t49 + t43 * t57 + (-t42 * t56 + t43 * t50) * t37) * t51 * t39) / (t51 ^ 2 * t39 + 0.1e1) 0, 0, 0, 0; 0 (t61 / t48 - t62 * t47 * t46) * t51 * t41, 0, 0, t74 * t41, 0;];
Ja_rot  = t1;
