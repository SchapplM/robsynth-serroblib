% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:27
% EndTime: 2019-02-26 19:46:27
% DurationCPUTime: 0.08s
% Computational Cost: add. (252->19), mult. (626->46), div. (35->9), fcn. (880->13), ass. (0->35)
t67 = sin(pkin(10));
t68 = sin(pkin(6));
t77 = t67 * t68;
t71 = cos(pkin(6));
t66 = sin(pkin(11));
t69 = cos(pkin(11));
t72 = sin(qJ(2));
t73 = cos(qJ(2));
t75 = t73 * t66 + t72 * t69;
t60 = t75 * t71;
t61 = t72 * t66 - t73 * t69;
t70 = cos(pkin(10));
t55 = -t67 * t60 - t70 * t61;
t65 = qJ(4) + pkin(12);
t63 = sin(t65);
t64 = cos(t65);
t50 = t55 * t64 + t63 * t77;
t48 = 0.1e1 / t50 ^ 2;
t49 = t55 * t63 - t64 * t77;
t76 = t49 ^ 2 * t48 + 0.1e1;
t74 = t61 * t71;
t59 = t75 * t68;
t58 = t61 * t68;
t57 = 0.1e1 / t58 ^ 2;
t53 = t67 * t74 - t70 * t75;
t52 = -t67 * t75 - t70 * t74;
t51 = -t70 * t60 + t67 * t61;
t47 = atan2(t52, t58);
t45 = cos(t47);
t44 = sin(t47);
t43 = 0.1e1 / t76;
t42 = t44 * t52 + t45 * t58;
t41 = 0.1e1 / t42 ^ 2;
t39 = (t51 / t58 - t59 * t52 * t57) / (t52 ^ 2 * t57 + 0.1e1);
t1 = [0, t39, 0, 0, 0, 0; 0 (t55 / t42 + (t44 * t51 + t45 * t59 + (-t44 * t58 + t45 * t52) * t39) * t53 * t41) / (t53 ^ 2 * t41 + 0.1e1) 0, 0, 0, 0; 0 (t63 / t50 - t64 * t49 * t48) * t53 * t43, 0, t76 * t43, 0, 0;];
Ja_rot  = t1;
