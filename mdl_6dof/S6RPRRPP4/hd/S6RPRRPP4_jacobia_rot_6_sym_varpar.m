% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:03
% EndTime: 2019-02-26 20:58:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (640->28), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->42)
t65 = pkin(9) + qJ(3);
t61 = sin(t65);
t82 = t61 ^ 2;
t63 = cos(t65);
t66 = qJ(4) + pkin(10);
t64 = cos(t66);
t69 = cos(qJ(1));
t71 = t69 * t64;
t62 = sin(t66);
t68 = sin(qJ(1));
t74 = t68 * t62;
t49 = t63 * t74 + t71;
t76 = t61 * t62;
t45 = atan2(-t49, t76);
t42 = sin(t45);
t43 = cos(t45);
t41 = -t42 * t49 + t43 * t76;
t40 = 0.1e1 / t41 ^ 2;
t72 = t69 * t62;
t73 = t68 * t64;
t52 = t63 * t72 - t73;
t81 = t40 * t52;
t79 = t43 * t49;
t78 = t52 ^ 2 * t40;
t57 = 0.1e1 / t61;
t59 = 0.1e1 / t62;
t77 = t57 * t59;
t75 = t61 * t69;
t53 = t63 * t71 + t74;
t48 = 0.1e1 / t53 ^ 2;
t70 = t69 ^ 2 * t82 * t48;
t60 = 0.1e1 / t62 ^ 2;
t58 = 0.1e1 / t82;
t51 = t63 * t73 - t72;
t47 = 0.1e1 / t53;
t46 = 0.1e1 / (0.1e1 + t70);
t44 = 0.1e1 / (t49 ^ 2 * t58 * t60 + 0.1e1);
t39 = 0.1e1 / t41;
t38 = (t49 * t58 * t59 * t63 + t68) * t44;
t37 = 0.1e1 / (0.1e1 + t78);
t36 = (t49 * t60 * t64 - t51 * t59) * t57 * t44;
t1 = [-t52 * t44 * t77, 0, t38, t36, 0, 0; (-t49 * t39 - (-t42 + (t77 * t79 + t42) * t44) * t78) * t37, 0 (t38 * t79 * t81 + (-t39 * t75 - (t43 * t63 + (-t38 + t68) * t61 * t42) * t81) * t62) * t37 (t53 * t39 - (t43 * t61 * t64 - t42 * t51 + (-t42 * t76 - t79) * t36) * t81) * t37, 0, 0; (-t48 * t51 * t69 + t47 * t68) * t61 * t46, 0 (-t47 * t63 * t69 - t64 * t70) * t46, -t52 * t48 * t46 * t75, 0, 0;];
Ja_rot  = t1;
