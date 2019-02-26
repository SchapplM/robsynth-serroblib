% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:03
% EndTime: 2019-02-26 21:09:03
% DurationCPUTime: 0.15s
% Computational Cost: add. (905->29), mult. (689->69), div. (139->11), fcn. (1046->9), ass. (0->42)
t79 = sin(qJ(3));
t91 = t79 ^ 2;
t74 = qJ(1) + pkin(10);
t70 = sin(t74);
t71 = cos(t74);
t78 = qJ(4) + qJ(5);
t73 = cos(t78);
t72 = sin(t78);
t80 = cos(qJ(3));
t84 = t72 * t80;
t60 = t70 * t84 + t71 * t73;
t82 = t79 * t72;
t57 = atan2(-t60, t82);
t53 = sin(t57);
t54 = cos(t57);
t51 = -t53 * t60 + t54 * t82;
t50 = 0.1e1 / t51 ^ 2;
t63 = -t70 * t73 + t71 * t84;
t90 = t50 * t63;
t89 = t50 * t63 ^ 2;
t87 = t54 * t60;
t68 = 0.1e1 / t72;
t76 = 0.1e1 / t79;
t86 = t68 * t76;
t85 = t71 * t79;
t83 = t73 * t80;
t64 = t70 * t72 + t71 * t83;
t59 = 0.1e1 / t64 ^ 2;
t81 = t59 * t71 ^ 2 * t91;
t77 = 0.1e1 / t91;
t69 = 0.1e1 / t72 ^ 2;
t62 = t70 * t83 - t71 * t72;
t58 = 0.1e1 / t64;
t56 = 0.1e1 / (t60 ^ 2 * t69 * t77 + 0.1e1);
t55 = 0.1e1 / (0.1e1 + t81);
t52 = t63 * t59 * t55 * t85;
t49 = 0.1e1 / t51;
t48 = (t60 * t68 * t77 * t80 + t70) * t56;
t47 = 0.1e1 / (0.1e1 + t89);
t46 = (t60 * t69 * t73 - t62 * t68) * t76 * t56;
t45 = (t64 * t49 - (t54 * t73 * t79 - t53 * t62 + (-t53 * t82 - t87) * t46) * t90) * t47;
t1 = [-t63 * t56 * t86, 0, t48, t46, t46, 0; (-t60 * t49 - (-t53 + (t86 * t87 + t53) * t56) * t89) * t47, 0 (t48 * t87 * t90 + (-t49 * t85 - (t54 * t80 + (-t48 + t70) * t79 * t53) * t90) * t72) * t47, t45, t45, 0; (-t59 * t62 * t71 + t58 * t70) * t79 * t55, 0 (-t58 * t71 * t80 - t73 * t81) * t55, -t52, -t52, 0;];
Ja_rot  = t1;
