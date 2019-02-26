% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:43
% EndTime: 2019-02-26 22:34:43
% DurationCPUTime: 0.10s
% Computational Cost: add. (339->26), mult. (503->68), div. (64->9), fcn. (719->11), ass. (0->42)
t88 = qJ(3) + qJ(4);
t84 = cos(t88);
t91 = sin(qJ(1));
t93 = cos(qJ(2));
t83 = sin(t88);
t94 = cos(qJ(1));
t97 = t94 * t83;
t76 = -t91 * t84 + t93 * t97;
t96 = t94 * t84;
t77 = t91 * t83 + t93 * t96;
t89 = sin(qJ(6));
t92 = cos(qJ(6));
t68 = t76 * t89 + t77 * t92;
t66 = 0.1e1 / t68 ^ 2;
t67 = -t76 * t92 + t77 * t89;
t105 = t66 * t67;
t104 = t67 ^ 2 * t66;
t90 = sin(qJ(2));
t99 = t91 * t90;
t81 = atan2(t99, t93);
t78 = sin(t81);
t79 = cos(t81);
t72 = t78 * t99 + t79 * t93;
t71 = 0.1e1 / t72 ^ 2;
t103 = t71 * t94 ^ 2;
t85 = t90 ^ 2;
t102 = t85 / t93 ^ 2;
t101 = t90 * t94;
t80 = 0.1e1 / (t91 ^ 2 * t102 + 0.1e1);
t100 = t91 * t80;
t98 = t91 * t93;
t95 = 0.1e1 + t104;
t86 = 0.1e1 / t93;
t75 = -t84 * t98 + t97;
t74 = -t83 * t98 - t96;
t73 = (0.1e1 + t102) * t100;
t70 = 0.1e1 / t72;
t69 = 0.1e1 / (t85 * t103 + 0.1e1);
t65 = 0.1e1 / t68;
t64 = 0.1e1 / t95;
t63 = (-t68 * t65 - t104) * t64;
t1 = [t86 * t80 * t101, t73, 0, 0, 0, 0; (t70 * t99 + (t79 * t85 * t86 * t100 + (-t80 + 0.1e1) * t90 * t78) * t90 * t103) * t69 (-t93 * t70 + (t78 * t98 - t79 * t90 + (-t78 * t93 + t79 * t99) * t73) * t90 * t71) * t94 * t69, 0, 0, 0, 0; ((-t74 * t92 + t75 * t89) * t65 - (t74 * t89 + t75 * t92) * t105) * t64 ((t83 * t92 - t84 * t89) * t65 - (-t83 * t89 - t84 * t92) * t105) * t64 * t101, t63, t63, 0, t95 * t64;];
Ja_rot  = t1;
