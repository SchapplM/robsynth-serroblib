% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR9_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobia_rot_3_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:23
% EndTime: 2019-02-26 20:53:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (156->21), mult. (450->48), div. (25->9), fcn. (637->13), ass. (0->38)
t61 = cos(pkin(6));
t59 = cos(pkin(12));
t65 = cos(qJ(1));
t69 = t65 * t59;
t56 = sin(pkin(12));
t63 = sin(qJ(1));
t72 = t63 * t56;
t51 = -t61 * t69 + t72;
t57 = sin(pkin(7));
t60 = cos(pkin(7));
t58 = sin(pkin(6));
t73 = t58 * t65;
t46 = -t51 * t57 + t60 * t73;
t50 = -t58 * t59 * t57 + t61 * t60;
t45 = atan2(t46, t50);
t42 = sin(t45);
t43 = cos(t45);
t37 = t42 * t46 + t43 * t50;
t70 = t65 * t56;
t71 = t63 * t59;
t53 = -t61 * t71 - t70;
t74 = t58 * t63;
t47 = t53 * t57 - t60 * t74;
t75 = t47 ^ 2 / t37 ^ 2;
t54 = -t61 * t72 + t69;
t62 = sin(qJ(3));
t64 = cos(qJ(3));
t66 = t53 * t60 + t57 * t74;
t41 = t54 * t64 + t66 * t62;
t39 = 0.1e1 / t41 ^ 2;
t40 = t54 * t62 - t66 * t64;
t68 = t40 ^ 2 * t39 + 0.1e1;
t67 = t51 * t60 + t57 * t73;
t52 = -t61 * t70 - t71;
t49 = 0.1e1 / t50;
t44 = 0.1e1 / (0.1e1 + t46 ^ 2 / t50 ^ 2);
t38 = 0.1e1 / t68;
t1 = [t47 * t49 * t44, 0, 0, 0, 0, 0; (t46 / t37 + (t42 + (t43 * t46 * t49 - t42) * t44) * t75) / (0.1e1 + t75) 0, 0, 0, 0, 0; ((t52 * t62 - t67 * t64) / t41 - (t52 * t64 + t67 * t62) * t40 * t39) * t38, 0, t68 * t38, 0, 0, 0;];
Ja_rot  = t1;
