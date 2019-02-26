% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:05
% EndTime: 2019-02-26 20:32:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (414->26), mult. (947->75), div. (100->11), fcn. (1427->11), ass. (0->42)
t68 = sin(qJ(4));
t87 = t68 ^ 2;
t86 = cos(qJ(1));
t85 = sin(qJ(1));
t74 = sin(pkin(9));
t75 = cos(pkin(9));
t56 = -t85 * t74 - t86 * t75;
t57 = t86 * t74 - t85 * t75;
t69 = cos(qJ(5));
t67 = sin(qJ(5));
t70 = cos(qJ(4));
t78 = t67 * t70;
t71 = -t56 * t69 + t57 * t78;
t77 = t68 * t67;
t44 = atan2(t71, -t77);
t42 = sin(t44);
t43 = cos(t44);
t40 = t42 * t71 - t43 * t77;
t39 = 0.1e1 / t40 ^ 2;
t51 = -t56 * t78 - t57 * t69;
t84 = t39 * t51;
t82 = t43 * t71;
t81 = t51 ^ 2 * t39;
t80 = t56 * t68;
t62 = 0.1e1 / t67;
t65 = 0.1e1 / t68;
t79 = t62 * t65;
t76 = t69 * t70;
t52 = -t56 * t76 + t57 * t67;
t47 = 0.1e1 / t52 ^ 2;
t73 = t56 ^ 2 * t87 * t47;
t72 = t56 * t67 + t57 * t76;
t66 = 0.1e1 / t87;
t63 = 0.1e1 / t67 ^ 2;
t46 = 0.1e1 / t52;
t45 = 0.1e1 / (t63 * t66 * t71 ^ 2 + 0.1e1);
t41 = 0.1e1 / (0.1e1 + t73);
t38 = 0.1e1 / t40;
t37 = (t62 * t66 * t70 * t71 + t57) * t45;
t36 = 0.1e1 / (0.1e1 + t81);
t35 = (t63 * t69 * t71 - t62 * t72) * t65 * t45;
t1 = [t51 * t45 * t79, 0, 0, t37, t35, 0; (t71 * t38 - (-t42 + (t79 * t82 + t42) * t45) * t81) * t36, 0, 0 (-t37 * t82 * t84 + (t38 * t80 - (-t43 * t70 + (t37 - t57) * t42 * t68) * t84) * t67) * t36 (t52 * t38 - (-t43 * t68 * t69 + t42 * t72 + (t42 * t77 + t82) * t35) * t84) * t36, 0; (-t72 * t56 * t47 - t57 * t46) * t68 * t41, 0, 0 (t46 * t56 * t70 - t69 * t73) * t41, t51 * t47 * t41 * t80, 0;];
Ja_rot  = t1;
