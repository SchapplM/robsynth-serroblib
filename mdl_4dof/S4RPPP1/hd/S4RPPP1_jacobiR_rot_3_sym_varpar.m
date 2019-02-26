% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S4RPPP1_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiR_rot_3_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiR_rot_3_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:29:48
% EndTime: 2019-02-26 19:29:48
% DurationCPUTime: 0.02s
% Computational Cost: add. (4->4), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
t45 = sin(pkin(6));
t49 = sin(qJ(1));
t54 = t49 * t45;
t47 = cos(pkin(6));
t53 = t49 * t47;
t50 = cos(qJ(1));
t52 = t50 * t45;
t51 = t50 * t47;
t48 = cos(pkin(4));
t46 = sin(pkin(4));
t1 = [t50 * t46, 0, 0, 0; t49 * t46, 0, 0, 0; 0, 0, 0, 0; t48 * t52 + t53, 0, 0, 0; t48 * t54 - t51, 0, 0, 0; 0, 0, 0, 0; t48 * t51 - t54, 0, 0, 0; t48 * t53 + t52, 0, 0, 0; 0, 0, 0, 0;];
JR_rot  = t1;
