% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:59
% EndTime: 2019-02-26 22:26:59
% DurationCPUTime: 0.11s
% Computational Cost: add. (243->21), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->35)
t61 = cos(qJ(2));
t59 = sin(qJ(2));
t60 = sin(qJ(1));
t67 = t60 * t59;
t51 = atan2(-t67, -t61);
t49 = sin(t51);
t50 = cos(t51);
t43 = -t49 * t67 - t50 * t61;
t42 = 0.1e1 / t43 ^ 2;
t62 = cos(qJ(1));
t72 = t42 * t62 ^ 2;
t55 = qJ(3) + qJ(4) + pkin(10);
t53 = sin(t55);
t54 = cos(t55);
t64 = t62 * t54;
t48 = t60 * t53 + t61 * t64;
t46 = 0.1e1 / t48 ^ 2;
t65 = t62 * t53;
t47 = -t60 * t54 + t61 * t65;
t71 = t46 * t47;
t56 = t59 ^ 2;
t70 = t56 / t61 ^ 2;
t69 = t59 * t62;
t52 = 0.1e1 / (t60 ^ 2 * t70 + 0.1e1);
t68 = t60 * t52;
t66 = t60 * t61;
t63 = t47 ^ 2 * t46 + 0.1e1;
t57 = 0.1e1 / t61;
t45 = 0.1e1 / t48;
t44 = (0.1e1 + t70) * t68;
t41 = 0.1e1 / t43;
t40 = 0.1e1 / (t56 * t72 + 0.1e1);
t39 = 0.1e1 / t63;
t38 = t63 * t39;
t1 = [t57 * t52 * t69, t44, 0, 0, 0, 0; (-t41 * t67 - (-t50 * t56 * t57 * t68 + (t52 - 0.1e1) * t59 * t49) * t59 * t72) * t40 (t61 * t41 - (-t49 * t66 + t50 * t59 + (t49 * t61 - t50 * t67) * t44) * t59 * t42) * t62 * t40, 0, 0, 0, 0; ((-t53 * t66 - t64) * t45 - (-t54 * t66 + t65) * t71) * t39 (-t45 * t53 + t54 * t71) * t39 * t69, t38, t38, 0, 0;];
Ja_rot  = t1;
