% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:38
% EndTime: 2019-02-26 20:10:38
% DurationCPUTime: 0.09s
% Computational Cost: add. (206->19), mult. (316->43), div. (52->11), fcn. (458->11), ass. (0->31)
t64 = sin(pkin(11));
t65 = sin(pkin(6));
t74 = t64 * t65;
t69 = cos(qJ(2));
t73 = t65 * t69;
t67 = cos(pkin(6));
t68 = sin(qJ(2));
t72 = t67 * t68;
t71 = t67 * t69;
t66 = cos(pkin(11));
t57 = -t64 * t72 + t66 * t69;
t61 = qJ(3) + qJ(4) + pkin(12);
t59 = sin(t61);
t60 = cos(t61);
t48 = t57 * t60 + t59 * t74;
t46 = 0.1e1 / t48 ^ 2;
t47 = t57 * t59 - t60 * t74;
t70 = t47 ^ 2 * t46 + 0.1e1;
t63 = 0.1e1 / t69 ^ 2;
t56 = t64 * t71 + t66 * t68;
t55 = t64 * t69 + t66 * t72;
t53 = t64 * t68 - t66 * t71;
t51 = atan2(-t53, -t73);
t50 = cos(t51);
t49 = sin(t51);
t45 = 0.1e1 / t70;
t44 = -t49 * t53 - t50 * t73;
t43 = 0.1e1 / t44 ^ 2;
t41 = (t55 / t69 + t68 * t53 * t63) / t65 / (0.1e1 + t53 ^ 2 / t65 ^ 2 * t63);
t40 = t70 * t45;
t1 = [0, t41, 0, 0, 0, 0; 0 (t57 / t44 - (t50 * t65 * t68 - t49 * t55 + (t49 * t73 - t50 * t53) * t41) * t56 * t43) / (t56 ^ 2 * t43 + 0.1e1) 0, 0, 0, 0; 0 (-t59 / t48 + t60 * t47 * t46) * t56 * t45, t40, t40, 0, 0;];
Ja_rot  = t1;
