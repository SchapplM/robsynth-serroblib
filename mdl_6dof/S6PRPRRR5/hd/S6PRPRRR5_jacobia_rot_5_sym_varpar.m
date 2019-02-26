% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:10
% EndTime: 2019-02-26 19:56:11
% DurationCPUTime: 0.10s
% Computational Cost: add. (150->18), mult. (316->43), div. (52->11), fcn. (458->11), ass. (0->31)
t58 = sin(pkin(11));
t59 = sin(pkin(6));
t68 = t58 * t59;
t62 = sin(qJ(2));
t67 = t59 * t62;
t61 = cos(pkin(6));
t66 = t61 * t62;
t63 = cos(qJ(2));
t65 = t61 * t63;
t60 = cos(pkin(11));
t51 = t58 * t65 + t60 * t62;
t57 = qJ(4) + qJ(5);
t53 = sin(t57);
t54 = cos(t57);
t43 = t51 * t53 + t54 * t68;
t41 = 0.1e1 / t43 ^ 2;
t42 = -t51 * t54 + t53 * t68;
t64 = t42 ^ 2 * t41 + 0.1e1;
t56 = 0.1e1 / t62 ^ 2;
t52 = -t58 * t66 + t60 * t63;
t49 = t58 * t63 + t60 * t66;
t48 = t58 * t62 - t60 * t65;
t47 = atan2(-t49, t67);
t45 = cos(t47);
t44 = sin(t47);
t40 = 0.1e1 / t64;
t39 = -t44 * t49 + t45 * t67;
t38 = 0.1e1 / t39 ^ 2;
t36 = (t48 / t62 + t63 * t49 * t56) / t59 / (0.1e1 + t49 ^ 2 / t59 ^ 2 * t56);
t35 = t64 * t40;
t1 = [0, t36, 0, 0, 0, 0; 0 (-t51 / t39 - (t45 * t59 * t63 + t44 * t48 + (-t44 * t67 - t45 * t49) * t36) * t52 * t38) / (t52 ^ 2 * t38 + 0.1e1) 0, 0, 0, 0; 0 (-t54 / t43 - t53 * t42 * t41) * t52 * t40, 0, t35, t35, 0;];
Ja_rot  = t1;
