% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-31 10:31
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiR_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiR_rot_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-31 10:31:13
% EndTime: 2019-05-31 10:31:13
% DurationCPUTime: 0.03s
% Computational Cost: add. (47->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t61 = qJ(2) + qJ(3);
t60 = cos(t61);
t62 = sin(qJ(4));
t72 = t60 * t62;
t63 = sin(qJ(1));
t71 = t63 * t62;
t64 = cos(qJ(4));
t70 = t63 * t64;
t65 = cos(qJ(1));
t69 = t65 * t62;
t68 = t65 * t64;
t59 = sin(t61);
t67 = t59 * t70;
t66 = t59 * t68;
t58 = t65 * t60;
t57 = t60 * t64;
t56 = t63 * t60;
t55 = t59 * t69;
t54 = t59 * t71;
t53 = t60 * t68 + t71;
t52 = -t60 * t69 + t70;
t51 = -t60 * t70 + t69;
t50 = t60 * t71 + t68;
t1 = [t51, -t66, -t66, t52, 0; t53, -t67, -t67, -t50, 0; 0, t57, t57, -t59 * t62, 0; t50, t55, t55, -t53, 0; t52, t54, t54, t51, 0; 0, -t72, -t72, -t59 * t64, 0; -t63 * t59, t58, t58, 0, 0; t65 * t59, t56, t56, 0, 0; 0, t59, t59, 0, 0;];
JR_rot  = t1;
