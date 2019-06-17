% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-12 14:37
% Revision: aab8d7cd0cba739f5e0ec8d53b8419901d1154b0 (2019-06-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiR_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiR_rot_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-12 14:37:34
% EndTime: 2019-06-12 14:37:34
% DurationCPUTime: 0.05s
% Computational Cost: add. (37->20), mult. (118->39), div. (0->0), fcn. (180->8), ass. (0->30)
t47 = sin(qJ(5));
t49 = sin(qJ(3));
t69 = t49 * t47;
t51 = cos(qJ(5));
t68 = t49 * t51;
t52 = cos(qJ(4));
t67 = t49 * t52;
t54 = cos(qJ(1));
t66 = t49 * t54;
t48 = sin(qJ(4));
t50 = sin(qJ(1));
t65 = t50 * t48;
t64 = t50 * t52;
t53 = cos(qJ(3));
t63 = t53 * t47;
t62 = t53 * t48;
t61 = t53 * t51;
t60 = t54 * t48;
t59 = t54 * t52;
t42 = t53 * t64 - t60;
t58 = -t42 * t47 + t50 * t68;
t57 = -t42 * t51 - t50 * t69;
t56 = -t51 * t67 + t63;
t55 = t47 * t67 + t61;
t44 = t53 * t59 + t65;
t43 = t53 * t60 - t64;
t41 = -t50 * t62 - t59;
t40 = t44 * t51 + t47 * t66;
t39 = -t44 * t47 + t51 * t66;
t1 = [t57, 0, t56 * t54, -t43 * t51, t39; t40, 0, t56 * t50, t41 * t51, t58; 0, 0, t52 * t61 + t69, -t48 * t68, -t55; -t58, 0, t55 * t54, t43 * t47, -t40; t39, 0, t55 * t50, -t41 * t47, t57; 0, 0, -t52 * t63 + t68, t48 * t69, t56; t41, 0, -t49 * t60, t44, 0; t43, 0, -t49 * t65, t42, 0; 0, 0, t62, t67, 0;];
JR_rot  = t1;
