% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:52:47
% EndTime: 2019-02-26 20:52:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (294->21), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
t63 = qJ(3) + pkin(10);
t59 = sin(t63);
t60 = cos(t63);
t66 = cos(qJ(1));
t70 = t66 * t60;
t54 = atan2(-t70, t59);
t52 = sin(t54);
t53 = cos(t54);
t45 = -t52 * t70 + t53 * t59;
t44 = 0.1e1 / t45 ^ 2;
t65 = sin(qJ(1));
t78 = t44 * t65 ^ 2;
t64 = qJ(5) + qJ(6);
t61 = sin(t64);
t69 = t66 * t61;
t62 = cos(t64);
t72 = t65 * t62;
t51 = t59 * t72 + t69;
t49 = 0.1e1 / t51 ^ 2;
t68 = t66 * t62;
t73 = t65 * t61;
t50 = t59 * t73 - t68;
t77 = t49 * t50;
t76 = t52 * t59;
t58 = t60 ^ 2;
t75 = 0.1e1 / t59 ^ 2 * t58;
t74 = t60 * t65;
t55 = 0.1e1 / (t66 ^ 2 * t75 + 0.1e1);
t71 = t66 * t55;
t67 = t50 ^ 2 * t49 + 0.1e1;
t56 = 0.1e1 / t59;
t48 = 0.1e1 / t51;
t47 = 0.1e1 / t67;
t46 = (0.1e1 + t75) * t71;
t43 = 0.1e1 / t45;
t42 = 0.1e1 / (t58 * t78 + 0.1e1);
t41 = t67 * t47;
t1 = [t56 * t55 * t74, 0, t46, 0, 0, 0; (-t43 * t70 + (-t53 * t56 * t58 * t71 + (-t55 + 0.1e1) * t60 * t52) * t60 * t78) * t42, 0 (t59 * t43 + (t66 * t76 + t53 * t60 + (-t53 * t70 - t76) * t46) * t60 * t44) * t65 * t42, 0, 0, 0; ((t59 * t69 + t72) * t48 - (t59 * t68 - t73) * t77) * t47, 0 (t48 * t61 - t62 * t77) * t47 * t74, 0, t41, t41;];
Ja_rot  = t1;
