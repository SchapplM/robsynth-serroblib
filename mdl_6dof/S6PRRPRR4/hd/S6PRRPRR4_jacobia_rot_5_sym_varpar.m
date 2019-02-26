% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:52
% EndTime: 2019-02-26 20:05:52
% DurationCPUTime: 0.10s
% Computational Cost: add. (183->21), mult. (524->53), div. (53->11), fcn. (745->13), ass. (0->35)
t55 = sin(pkin(11));
t57 = cos(pkin(11));
t64 = cos(qJ(2));
t58 = cos(pkin(6));
t61 = sin(qJ(2));
t67 = t58 * t61;
t50 = -t55 * t67 + t57 * t64;
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t56 = sin(pkin(6));
t69 = t55 * t56;
t41 = t50 * t60 - t63 * t69;
t42 = t50 * t63 + t60 * t69;
t59 = sin(qJ(5));
t62 = cos(qJ(5));
t40 = t41 * t59 + t42 * t62;
t38 = 0.1e1 / t40 ^ 2;
t39 = -t41 * t62 + t42 * t59;
t70 = t39 ^ 2 * t38;
t68 = t56 * t64;
t66 = t58 * t64;
t65 = 0.1e1 + t70;
t54 = 0.1e1 / t64 ^ 2;
t49 = -t55 * t66 - t57 * t61;
t48 = t55 * t64 + t57 * t67;
t47 = t55 * t61 - t57 * t66;
t46 = atan2(t47, t68);
t44 = cos(t46);
t43 = sin(t46);
t37 = 0.1e1 / t40;
t36 = t43 * t47 + t44 * t68;
t35 = 0.1e1 / t36 ^ 2;
t33 = (t48 / t64 + t61 * t47 * t54) / t56 / (0.1e1 + t47 ^ 2 / t56 ^ 2 * t54);
t32 = 0.1e1 / t65;
t1 = [0, t33, 0, 0, 0, 0; 0 (-t50 / t36 - (-t44 * t56 * t61 + t43 * t48 + (-t43 * t68 + t44 * t47) * t33) * t49 * t35) / (t49 ^ 2 * t35 + 0.1e1) 0, 0, 0, 0; 0 ((t59 * t63 - t60 * t62) * t37 - (t59 * t60 + t62 * t63) * t39 * t38) * t32 * t49 (-t40 * t37 - t70) * t32, 0, t65 * t32, 0;];
Ja_rot  = t1;
