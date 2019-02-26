% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR12
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR12_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:02
% EndTime: 2019-02-26 20:55:02
% DurationCPUTime: 0.11s
% Computational Cost: add. (136->20), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->36)
t57 = cos(qJ(3));
t55 = sin(qJ(3));
t58 = cos(qJ(1));
t60 = t58 * t55;
t48 = atan2(t60, t57);
t45 = sin(t48);
t46 = cos(t48);
t39 = t45 * t60 + t46 * t57;
t38 = 0.1e1 / t39 ^ 2;
t56 = sin(qJ(1));
t69 = t38 * t56 ^ 2;
t54 = qJ(5) + qJ(6);
t49 = sin(t54);
t50 = cos(t54);
t61 = t58 * t50;
t64 = t56 * t57;
t44 = -t49 * t64 + t61;
t42 = 0.1e1 / t44 ^ 2;
t62 = t58 * t49;
t43 = t50 * t64 + t62;
t68 = t42 * t43;
t67 = t45 * t57;
t51 = t55 ^ 2;
t66 = t51 / t57 ^ 2;
t65 = t55 * t56;
t47 = 0.1e1 / (t58 ^ 2 * t66 + 0.1e1);
t63 = t58 * t47;
t59 = t43 ^ 2 * t42 + 0.1e1;
t52 = 0.1e1 / t57;
t41 = 0.1e1 / t44;
t40 = (0.1e1 + t66) * t63;
t37 = 0.1e1 / t39;
t36 = 0.1e1 / t59;
t35 = 0.1e1 / (t51 * t69 + 0.1e1);
t34 = t59 * t36;
t1 = [-t52 * t47 * t65, 0, t40, 0, 0, 0; (t37 * t60 - (-t46 * t51 * t52 * t63 + (t47 - 0.1e1) * t55 * t45) * t55 * t69) * t35, 0 (t57 * t37 - (t58 * t67 - t46 * t55 + (t46 * t60 - t67) * t40) * t55 * t38) * t56 * t35, 0, 0, 0; ((-t56 * t49 + t57 * t61) * t41 - (-t56 * t50 - t57 * t62) * t68) * t36, 0 (-t41 * t50 - t49 * t68) * t36 * t65, 0, t34, t34;];
Ja_rot  = t1;
