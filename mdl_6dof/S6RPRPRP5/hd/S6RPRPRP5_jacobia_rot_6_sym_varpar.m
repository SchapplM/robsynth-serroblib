% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP5
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
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:52
% EndTime: 2019-02-26 20:45:53
% DurationCPUTime: 0.13s
% Computational Cost: add. (640->28), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->42)
t62 = pkin(9) + qJ(3);
t58 = sin(t62);
t78 = t58 ^ 2;
t60 = cos(t62);
t61 = pkin(10) + qJ(5);
t59 = cos(t61);
t65 = cos(qJ(1));
t67 = t65 * t59;
t57 = sin(t61);
t64 = sin(qJ(1));
t70 = t64 * t57;
t45 = t60 * t70 + t67;
t72 = t58 * t57;
t41 = atan2(-t45, t72);
t38 = sin(t41);
t39 = cos(t41);
t37 = -t38 * t45 + t39 * t72;
t36 = 0.1e1 / t37 ^ 2;
t68 = t65 * t57;
t69 = t64 * t59;
t48 = t60 * t68 - t69;
t77 = t36 * t48;
t75 = t39 * t45;
t74 = t48 ^ 2 * t36;
t52 = 0.1e1 / t57;
t55 = 0.1e1 / t58;
t73 = t52 * t55;
t71 = t58 * t65;
t49 = t60 * t67 + t70;
t44 = 0.1e1 / t49 ^ 2;
t66 = t65 ^ 2 * t78 * t44;
t56 = 0.1e1 / t78;
t53 = 0.1e1 / t57 ^ 2;
t47 = t60 * t69 - t68;
t43 = 0.1e1 / t49;
t42 = 0.1e1 / (0.1e1 + t66);
t40 = 0.1e1 / (t45 ^ 2 * t56 * t53 + 0.1e1);
t35 = 0.1e1 / t37;
t34 = (t45 * t52 * t56 * t60 + t64) * t40;
t33 = 0.1e1 / (0.1e1 + t74);
t32 = (t45 * t53 * t59 - t47 * t52) * t55 * t40;
t1 = [-t48 * t40 * t73, 0, t34, 0, t32, 0; (-t45 * t35 - (-t38 + (t73 * t75 + t38) * t40) * t74) * t33, 0 (t34 * t75 * t77 + (-t35 * t71 - (t39 * t60 + (-t34 + t64) * t58 * t38) * t77) * t57) * t33, 0 (t49 * t35 - (t39 * t58 * t59 - t38 * t47 + (-t38 * t72 - t75) * t32) * t77) * t33, 0; (-t44 * t47 * t65 + t43 * t64) * t58 * t42, 0 (-t43 * t60 * t65 - t59 * t66) * t42, 0, -t48 * t44 * t42 * t71, 0;];
Ja_rot  = t1;
