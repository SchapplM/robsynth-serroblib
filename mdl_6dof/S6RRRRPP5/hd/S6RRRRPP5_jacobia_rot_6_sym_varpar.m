% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:27:38
% EndTime: 2019-02-26 22:27:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (141->20), mult. (257->55), div. (59->9), fcn. (377->9), ass. (0->35)
t70 = cos(qJ(2));
t68 = sin(qJ(2));
t69 = sin(qJ(1));
t75 = t69 * t68;
t61 = atan2(t75, t70);
t58 = sin(t61);
t59 = cos(t61);
t51 = t58 * t75 + t59 * t70;
t50 = 0.1e1 / t51 ^ 2;
t71 = cos(qJ(1));
t81 = t50 * t71 ^ 2;
t67 = qJ(3) + qJ(4);
t62 = sin(t67);
t63 = cos(t67);
t72 = t71 * t63;
t57 = t69 * t62 + t70 * t72;
t55 = 0.1e1 / t57 ^ 2;
t73 = t71 * t62;
t56 = t69 * t63 - t70 * t73;
t80 = t56 ^ 2 * t55;
t79 = t55 * t56;
t64 = t68 ^ 2;
t78 = t64 / t70 ^ 2;
t77 = t68 * t71;
t60 = 0.1e1 / (t69 ^ 2 * t78 + 0.1e1);
t76 = t69 * t60;
t74 = t69 * t70;
t65 = 0.1e1 / t70;
t54 = 0.1e1 / t57;
t52 = (0.1e1 + t78) * t76;
t49 = 0.1e1 / t51;
t48 = 0.1e1 / (0.1e1 + t80);
t47 = 0.1e1 / (t64 * t81 + 0.1e1);
t46 = (-t54 * t57 - t80) * t48;
t1 = [t65 * t60 * t77, t52, 0, 0, 0, 0; (t49 * t75 + (t59 * t64 * t65 * t76 + (-t60 + 0.1e1) * t68 * t58) * t68 * t81) * t47 (-t70 * t49 + (t58 * t74 - t59 * t68 + (-t58 * t70 + t59 * t75) * t52) * t68 * t50) * t71 * t47, 0, 0, 0, 0; ((t62 * t74 + t72) * t54 - (-t63 * t74 + t73) * t79) * t48 (t54 * t62 + t63 * t79) * t48 * t77, t46, t46, 0, 0;];
Ja_rot  = t1;
