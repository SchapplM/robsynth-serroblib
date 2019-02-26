% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:15
% EndTime: 2019-02-26 20:44:15
% DurationCPUTime: 0.13s
% Computational Cost: add. (572->28), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->42)
t60 = qJ(3) + pkin(10);
t56 = sin(t60);
t78 = t56 ^ 2;
t58 = cos(t60);
t61 = qJ(1) + pkin(9);
t59 = cos(t61);
t65 = cos(qJ(5));
t67 = t59 * t65;
t57 = sin(t61);
t64 = sin(qJ(5));
t70 = t57 * t64;
t45 = t58 * t70 + t67;
t71 = t56 * t64;
t42 = atan2(-t45, t71);
t39 = sin(t42);
t40 = cos(t42);
t37 = -t39 * t45 + t40 * t71;
t36 = 0.1e1 / t37 ^ 2;
t68 = t59 * t64;
t69 = t57 * t65;
t48 = t58 * t68 - t69;
t77 = t36 * t48;
t75 = t40 * t45;
t74 = t48 ^ 2 * t36;
t53 = 0.1e1 / t56;
t62 = 0.1e1 / t64;
t73 = t53 * t62;
t72 = t56 * t59;
t49 = t58 * t67 + t70;
t44 = 0.1e1 / t49 ^ 2;
t66 = t59 ^ 2 * t78 * t44;
t63 = 0.1e1 / t64 ^ 2;
t54 = 0.1e1 / t78;
t47 = t58 * t69 - t68;
t43 = 0.1e1 / t49;
t41 = 0.1e1 / (t45 ^ 2 * t54 * t63 + 0.1e1);
t38 = 0.1e1 / (0.1e1 + t66);
t35 = 0.1e1 / t37;
t34 = (t45 * t54 * t58 * t62 + t57) * t41;
t33 = 0.1e1 / (0.1e1 + t74);
t32 = (t45 * t63 * t65 - t47 * t62) * t53 * t41;
t1 = [-t48 * t41 * t73, 0, t34, 0, t32, 0; (-t45 * t35 - (-t39 + (t73 * t75 + t39) * t41) * t74) * t33, 0 (t34 * t75 * t77 + (-t35 * t72 - (t40 * t58 + (-t34 + t57) * t56 * t39) * t77) * t64) * t33, 0 (t49 * t35 - (t40 * t56 * t65 - t39 * t47 + (-t39 * t71 - t75) * t32) * t77) * t33, 0; (-t44 * t47 * t59 + t43 * t57) * t56 * t38, 0 (-t43 * t58 * t59 - t65 * t66) * t38, 0, -t48 * t44 * t38 * t72, 0;];
Ja_rot  = t1;
