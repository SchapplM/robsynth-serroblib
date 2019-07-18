% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRR1_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiR_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiR_rot_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:25
% EndTime: 2019-07-18 13:26:25
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t57 = sin(qJ(3));
t58 = sin(qJ(1));
t68 = t58 * t57;
t59 = cos(qJ(4));
t67 = t58 * t59;
t56 = sin(qJ(4));
t60 = cos(qJ(3));
t66 = t60 * t56;
t65 = t60 * t59;
t61 = cos(qJ(1));
t64 = t61 * t57;
t63 = t61 * t59;
t62 = t61 * t60;
t55 = t58 * t56 + t59 * t62;
t54 = -t56 * t62 + t67;
t53 = t61 * t56 - t58 * t65;
t52 = t58 * t66 + t63;
t1 = [t53, 0, -t57 * t63, t54, 0; t55, 0, -t57 * t67, -t52, 0; 0, 0, t65, -t57 * t56, 0; t52, 0, t56 * t64, -t55, 0; t54, 0, t56 * t68, t53, 0; 0, 0, -t66, -t57 * t59, 0; -t68, 0, t62, 0, 0; t64, 0, t58 * t60, 0, 0; 0, 0, t57, 0, 0;];
JR_rot  = t1;
