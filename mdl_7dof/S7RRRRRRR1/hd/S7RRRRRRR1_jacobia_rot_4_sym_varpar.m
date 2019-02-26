% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Ja_rot [3x7]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S7RRRRRRR1_jacobia_rot_4_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_rot_4_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_rot_4_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:53
% EndTime: 2018-11-26 21:20:53
% DurationCPUTime: 0.14s
% Computational Cost: add. (252->32), mult. (662->85), div. (108->11), fcn. (985->11), ass. (0->45)
t56 = sin(qJ(3));
t60 = cos(qJ(3));
t62 = cos(qJ(1));
t64 = t62 * t60;
t58 = sin(qJ(1));
t61 = cos(qJ(2));
t66 = t58 * t61;
t45 = t56 * t66 - t64;
t57 = sin(qJ(2));
t70 = t57 * t56;
t43 = atan2(t45, -t70);
t41 = sin(t43);
t42 = cos(t43);
t35 = t41 * t45 - t42 * t70;
t33 = 0.1e1 / t35 ^ 2;
t65 = t62 * t56;
t47 = t58 * t60 + t61 * t65;
t75 = t33 * t47;
t49 = -t58 * t56 + t61 * t64;
t55 = sin(qJ(4));
t59 = cos(qJ(4));
t67 = t57 * t62;
t40 = t49 * t59 + t55 * t67;
t38 = 0.1e1 / t40 ^ 2;
t39 = t49 * t55 - t59 * t67;
t74 = t38 * t39;
t73 = t42 * t45;
t72 = t47 ^ 2 * t33;
t51 = 0.1e1 / t56;
t53 = 0.1e1 / t57;
t71 = t51 * t53;
t69 = t57 * t58;
t68 = t57 * t60;
t63 = t39 ^ 2 * t38 + 0.1e1;
t54 = 0.1e1 / t57 ^ 2;
t52 = 0.1e1 / t56 ^ 2;
t46 = t60 * t66 + t65;
t44 = 0.1e1 / (t45 ^ 2 * t54 * t52 + 0.1e1);
t37 = 0.1e1 / t40;
t36 = 0.1e1 / t63;
t34 = (t45 * t51 * t54 * t61 + t58) * t44;
t32 = 0.1e1 / t35;
t31 = 0.1e1 / (0.1e1 + t72);
t30 = (t45 * t52 * t60 - t46 * t51) * t53 * t44;
t1 = [-t47 * t44 * t71, t34, t30, 0, 0, 0, 0; (t45 * t32 + (t41 + (-t71 * t73 - t41) * t44) * t72) * t31 (t34 * t73 * t75 + (t32 * t67 + (-t42 * t61 + (t34 * t57 - t69) * t41) * t75) * t56) * t31 (-t49 * t32 + (-t42 * t68 + t41 * t46 + (t41 * t70 + t73) * t30) * t75) * t31, 0, 0, 0, 0; ((-t46 * t55 + t59 * t69) * t37 - (-t46 * t59 - t55 * t69) * t74) * t36 ((-t55 * t68 - t59 * t61) * t37 - (t55 * t61 - t59 * t68) * t74) * t36 * t62 (-t55 * t37 + t59 * t74) * t47 * t36, t63 * t36, 0, 0, 0;];
Ja_rot  = t1;
