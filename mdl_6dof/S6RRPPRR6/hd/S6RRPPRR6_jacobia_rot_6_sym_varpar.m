% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:25
% EndTime: 2019-02-26 21:31:25
% DurationCPUTime: 0.17s
% Computational Cost: add. (737->24), mult. (876->55), div. (85->9), fcn. (1282->11), ass. (0->38)
t85 = pkin(10) + qJ(5);
t71 = sin(t85);
t75 = cos(qJ(2));
t82 = cos(t85);
t93 = sin(qJ(2));
t65 = -t75 * t71 + t93 * t82;
t73 = sin(qJ(1));
t57 = t65 * t73;
t64 = t93 * t71 + t75 * t82;
t52 = atan2(t57, t64);
t49 = sin(t52);
t50 = cos(t52);
t47 = t49 * t57 + t50 * t64;
t46 = 0.1e1 / t47 ^ 2;
t94 = cos(qJ(1));
t61 = t65 * t94;
t87 = t61 ^ 2 * t46;
t44 = 0.1e1 / (0.1e1 + t87);
t45 = 0.1e1 / t47;
t59 = t64 * t73;
t60 = t64 * t94;
t90 = t50 * t57;
t63 = 0.1e1 / t64 ^ 2;
t51 = 0.1e1 / (t57 ^ 2 * t63 + 0.1e1);
t62 = 0.1e1 / t64;
t96 = (-t57 * t63 * t65 - t59 * t62) * t51;
t99 = ((-t49 * t59 + t50 * t65 + (-t49 * t64 + t90) * t96) * t46 * t61 + t60 * t45) * t44;
t72 = sin(qJ(6));
t74 = cos(qJ(6));
t56 = t60 * t74 - t73 * t72;
t54 = 0.1e1 / t56 ^ 2;
t55 = t60 * t72 + t73 * t74;
t84 = t55 ^ 2 * t54 + 0.1e1;
t48 = 0.1e1 / t84;
t53 = 0.1e1 / t56;
t89 = t54 * t55;
t95 = (-t72 * t53 + t74 * t89) * t48 * t61;
t1 = [t61 * t62 * t51, -t96, 0, 0, t96, 0; (t57 * t45 - (-t49 + (-t62 * t90 + t49) * t51) * t87) * t44, -t99, 0, 0, t99, 0; ((-t59 * t72 + t94 * t74) * t53 - (-t59 * t74 - t94 * t72) * t89) * t48, t95, 0, 0, -t95, t84 * t48;];
Ja_rot  = t1;
