% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:25
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14_jacobia_rot_3_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_rot_3_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_rot_3_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-03 10:25:33
% EndTime: 2019-01-03 10:25:33
% DurationCPUTime: 0.14s
% Computational Cost: add. (273->29), mult. (795->75), div. (50->9), fcn. (1135->13), ass. (0->47)
t61 = cos(pkin(6));
t64 = cos(qJ(2));
t65 = cos(qJ(1));
t68 = t64 * t65;
t62 = sin(qJ(2));
t63 = sin(qJ(1));
t70 = t63 * t62;
t50 = -t61 * t68 + t70;
t57 = sin(pkin(7));
t60 = cos(pkin(7));
t58 = sin(pkin(6));
t72 = t58 * t65;
t44 = -t50 * t57 + t60 * t72;
t49 = -t57 * t58 * t64 + t60 * t61;
t43 = atan2(t44, t49);
t40 = sin(t43);
t41 = cos(t43);
t34 = t40 * t44 + t41 * t49;
t33 = 0.1e1 / t34 ^ 2;
t69 = t63 * t64;
t71 = t62 * t65;
t52 = -t61 * t69 - t71;
t73 = t58 * t63;
t45 = t52 * t57 - t60 * t73;
t78 = t33 * t45 ^ 2;
t53 = -t61 * t70 + t68;
t56 = sin(pkin(14));
t59 = cos(pkin(14));
t66 = t52 * t60 + t57 * t73;
t39 = t53 * t59 + t66 * t56;
t37 = 0.1e1 / t39 ^ 2;
t38 = t53 * t56 - t66 * t59;
t77 = t37 * t38;
t76 = t41 * t44;
t75 = t53 * t60;
t74 = t58 * t62;
t67 = t50 * t60 + t57 * t72;
t51 = -t61 * t71 - t69;
t48 = 0.1e1 / t49 ^ 2;
t47 = 0.1e1 / t49;
t42 = 0.1e1 / (t44 ^ 2 * t48 + 0.1e1);
t36 = 0.1e1 / t39;
t35 = 0.1e1 / (t37 * t38 ^ 2 + 0.1e1);
t32 = 0.1e1 / t34;
t31 = 0.1e1 / (0.1e1 + t78);
t30 = (-t44 * t48 * t74 + t47 * t51) * t57 * t42;
t1 = [t45 * t47 * t42, t30, 0, 0, 0, 0; (t44 * t32 + (t40 + (t47 * t76 - t40) * t42) * t78) * t31 (t53 * t57 * t32 + ((t40 * t51 + t41 * t74) * t57 + (-t40 * t49 + t76) * t30) * t45 * t33) * t31, 0, 0, 0, 0; ((t51 * t56 - t67 * t59) * t36 - (t51 * t59 + t67 * t56) * t77) * t35 ((t52 * t56 + t59 * t75) * t36 - (t52 * t59 - t56 * t75) * t77) * t35, 0, 0, 0, 0;];
Ja_rot  = t1;
