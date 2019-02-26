% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR8_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:40
% EndTime: 2019-02-26 21:18:40
% DurationCPUTime: 0.11s
% Computational Cost: add. (323->20), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t63 = qJ(3) + qJ(4);
t61 = sin(t63);
t62 = cos(t63);
t67 = cos(qJ(1));
t71 = t67 * t62;
t56 = atan2(-t71, t61);
t54 = sin(t56);
t55 = cos(t56);
t47 = -t54 * t71 + t55 * t61;
t46 = 0.1e1 / t47 ^ 2;
t65 = sin(qJ(1));
t79 = t46 * t65 ^ 2;
t64 = sin(qJ(5));
t70 = t67 * t64;
t66 = cos(qJ(5));
t73 = t65 * t66;
t53 = t61 * t73 + t70;
t51 = 0.1e1 / t53 ^ 2;
t69 = t67 * t66;
t74 = t65 * t64;
t52 = t61 * t74 - t69;
t78 = t51 * t52;
t77 = t54 * t61;
t60 = t62 ^ 2;
t76 = 0.1e1 / t61 ^ 2 * t60;
t75 = t62 * t65;
t57 = 0.1e1 / (t67 ^ 2 * t76 + 0.1e1);
t72 = t67 * t57;
t68 = t52 ^ 2 * t51 + 0.1e1;
t58 = 0.1e1 / t61;
t50 = 0.1e1 / t53;
t49 = 0.1e1 / t68;
t48 = (0.1e1 + t76) * t72;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (t60 * t79 + 0.1e1);
t43 = (t50 * t64 - t66 * t78) * t49 * t75;
t42 = (t61 * t45 + (t67 * t77 + t55 * t62 + (-t55 * t71 - t77) * t48) * t62 * t46) * t65 * t44;
t1 = [t58 * t57 * t75, 0, t48, t48, 0, 0; (-t45 * t71 + (-t55 * t58 * t60 * t72 + (-t57 + 0.1e1) * t62 * t54) * t62 * t79) * t44, 0, t42, t42, 0, 0; ((t61 * t70 + t73) * t50 - (t61 * t69 - t74) * t78) * t49, 0, t43, t43, t68 * t49, 0;];
Ja_rot  = t1;
