% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR4_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:10
% EndTime: 2019-02-26 21:30:10
% DurationCPUTime: 0.13s
% Computational Cost: add. (289->21), mult. (821->56), div. (53->11), fcn. (1181->11), ass. (0->37)
t64 = sin(qJ(1));
t66 = cos(qJ(1));
t60 = sin(pkin(11));
t62 = cos(pkin(11));
t65 = cos(qJ(2));
t71 = cos(pkin(6));
t69 = t65 * t71;
t63 = sin(qJ(2));
t70 = t63 * t71;
t67 = -t60 * t70 + t62 * t69;
t68 = t65 * t60 + t63 * t62;
t44 = -t64 * t68 + t66 * t67;
t54 = t63 * t60 - t65 * t62;
t61 = sin(pkin(6));
t51 = t54 * t61;
t41 = atan2(t44, t51);
t39 = cos(t41);
t74 = t39 * t44;
t53 = t60 * t69 + t62 * t70;
t47 = -t64 * t53 - t66 * t54;
t59 = 0.1e1 / t64 ^ 2;
t73 = 0.1e1 / (0.1e1 + t47 ^ 2 * t59 / t61 ^ 2) / t61;
t38 = sin(t41);
t37 = t38 * t44 + t39 * t51;
t36 = 0.1e1 / t37 ^ 2;
t45 = -t64 * t67 - t66 * t68;
t72 = t45 ^ 2 * t36;
t43 = -t66 * t53 + t64 * t54;
t58 = 0.1e1 / t64;
t52 = t68 * t61;
t50 = 0.1e1 / t51 ^ 2;
t49 = 0.1e1 / t51;
t40 = 0.1e1 / (t44 ^ 2 * t50 + 0.1e1);
t35 = 0.1e1 / t37;
t34 = 0.1e1 / (0.1e1 + t72);
t33 = (-t44 * t50 * t52 + t43 * t49) * t40;
t1 = [t45 * t49 * t40, t33, 0, 0, 0, 0; (t44 * t35 + (t38 + (t49 * t74 - t38) * t40) * t72) * t34 (t47 * t35 + (t38 * t43 + t39 * t52 + (-t38 * t51 + t74) * t33) * t45 * t36) * t34, 0, 0, 0, 0; (-t66 * t47 * t59 + t43 * t58) * t73, t45 * t58 * t73, 0, 0, 0, 0;];
Ja_rot  = t1;
