% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:36
% EndTime: 2019-02-26 20:18:36
% DurationCPUTime: 0.09s
% Computational Cost: add. (251->19), mult. (364->43), div. (57->11), fcn. (521->11), ass. (0->31)
t63 = sin(pkin(12));
t64 = sin(pkin(6));
t73 = t63 * t64;
t68 = cos(qJ(2));
t72 = t64 * t68;
t66 = cos(pkin(6));
t67 = sin(qJ(2));
t71 = t66 * t67;
t70 = t66 * t68;
t65 = cos(pkin(12));
t56 = -t63 * t71 + t65 * t68;
t60 = qJ(3) + qJ(4) + qJ(5);
t58 = sin(t60);
t59 = cos(t60);
t47 = t56 * t59 + t58 * t73;
t45 = 0.1e1 / t47 ^ 2;
t46 = t56 * t58 - t59 * t73;
t69 = t46 ^ 2 * t45 + 0.1e1;
t62 = 0.1e1 / t68 ^ 2;
t55 = t63 * t70 + t65 * t67;
t54 = t63 * t68 + t65 * t71;
t52 = t63 * t67 - t65 * t70;
t50 = atan2(-t52, -t72);
t49 = cos(t50);
t48 = sin(t50);
t44 = 0.1e1 / t69;
t43 = -t48 * t52 - t49 * t72;
t42 = 0.1e1 / t43 ^ 2;
t40 = (t54 / t68 + t67 * t52 * t62) / t64 / (0.1e1 + t52 ^ 2 / t64 ^ 2 * t62);
t39 = t69 * t44;
t1 = [0, t40, 0, 0, 0, 0; 0 (t56 / t43 - (t49 * t64 * t67 - t48 * t54 + (t48 * t72 - t49 * t52) * t40) * t55 * t42) / (t55 ^ 2 * t42 + 0.1e1) 0, 0, 0, 0; 0 (-t58 / t47 + t59 * t46 * t45) * t55 * t44, t39, t39, t39, 0;];
Ja_rot  = t1;
