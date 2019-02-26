% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:15
% EndTime: 2019-02-26 20:51:15
% DurationCPUTime: 0.17s
% Computational Cost: add. (737->24), mult. (876->55), div. (85->9), fcn. (1282->11), ass. (0->38)
t84 = pkin(10) + qJ(3);
t70 = cos(t84);
t72 = sin(qJ(5));
t81 = sin(t84);
t92 = cos(qJ(5));
t64 = -t70 * t72 + t81 * t92;
t73 = sin(qJ(1));
t56 = t64 * t73;
t63 = t70 * t92 + t81 * t72;
t51 = atan2(t56, t63);
t48 = sin(t51);
t49 = cos(t51);
t46 = t48 * t56 + t49 * t63;
t45 = 0.1e1 / t46 ^ 2;
t93 = cos(qJ(1));
t60 = t64 * t93;
t86 = t60 ^ 2 * t45;
t43 = 0.1e1 / (0.1e1 + t86);
t44 = 0.1e1 / t46;
t58 = t63 * t73;
t59 = t63 * t93;
t89 = t49 * t56;
t62 = 0.1e1 / t63 ^ 2;
t50 = 0.1e1 / (t56 ^ 2 * t62 + 0.1e1);
t61 = 0.1e1 / t63;
t95 = (-t56 * t62 * t64 - t58 * t61) * t50;
t98 = (t45 * t60 * (-t48 * t58 + t49 * t64 + (-t48 * t63 + t89) * t95) + t59 * t44) * t43;
t71 = sin(qJ(6));
t74 = cos(qJ(6));
t55 = t59 * t74 - t73 * t71;
t53 = 0.1e1 / t55 ^ 2;
t54 = t59 * t71 + t73 * t74;
t83 = t54 ^ 2 * t53 + 0.1e1;
t47 = 0.1e1 / t83;
t52 = 0.1e1 / t55;
t88 = t53 * t54;
t94 = (-t71 * t52 + t74 * t88) * t47 * t60;
t1 = [t60 * t61 * t50, 0, -t95, 0, t95, 0; (t56 * t44 - (-t48 + (-t61 * t89 + t48) * t50) * t86) * t43, 0, -t98, 0, t98, 0; ((-t58 * t71 + t93 * t74) * t52 - (-t58 * t74 - t93 * t71) * t88) * t47, 0, t94, 0, -t94, t83 * t47;];
Ja_rot  = t1;
