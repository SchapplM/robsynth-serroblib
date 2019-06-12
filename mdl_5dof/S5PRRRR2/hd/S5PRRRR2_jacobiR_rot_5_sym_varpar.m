% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5PRRRR2
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
%   pkin=[a4]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-03 15:11
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_jacobiR_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5PRRRR2_jacobiR_rot_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-03 15:11:34
% EndTime: 2019-06-03 15:11:34
% DurationCPUTime: 0.05s
% Computational Cost: add. (47->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t51 = qJ(3) + qJ(4);
t50 = cos(t51);
t54 = cos(qJ(5));
t62 = t50 * t54;
t52 = sin(qJ(5));
t53 = sin(qJ(2));
t61 = t53 * t52;
t60 = t53 * t54;
t55 = cos(qJ(2));
t59 = t55 * t52;
t58 = t55 * t54;
t49 = sin(t51);
t57 = t49 * t60;
t56 = t49 * t58;
t48 = t55 * t50;
t47 = t50 * t52;
t46 = t53 * t50;
t45 = t49 * t59;
t44 = t49 * t61;
t43 = t50 * t58 + t61;
t42 = -t50 * t59 + t60;
t41 = -t50 * t60 + t59;
t40 = t50 * t61 + t58;
t1 = [0, t41, -t56, -t56, t42; 0, 0, -t62, -t62, t49 * t52; 0, t43, -t57, -t57, -t40; 0, t40, t45, t45, -t43; 0, 0, t47, t47, t49 * t54; 0, t42, t44, t44, t41; 0, -t53 * t49, t48, t48, 0; 0, 0, -t49, -t49, 0; 0, t55 * t49, t46, t46, 0;];
JR_rot  = t1;
