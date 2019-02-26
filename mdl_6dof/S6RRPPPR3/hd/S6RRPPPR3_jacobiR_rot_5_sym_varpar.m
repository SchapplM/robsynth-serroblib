% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:03
% EndTime: 2019-02-26 21:23:03
% DurationCPUTime: 0.03s
% Computational Cost: add. (10->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->11)
t49 = sin(qJ(2));
t50 = sin(qJ(1));
t56 = t50 * t49;
t51 = cos(qJ(2));
t55 = t50 * t51;
t52 = cos(qJ(1));
t54 = t52 * t49;
t53 = t52 * t51;
t48 = cos(pkin(9));
t47 = sin(pkin(9));
t1 = [-t52 * t47 - t48 * t56, t48 * t53, 0, 0, 0, 0; -t50 * t47 + t48 * t54, t48 * t55, 0, 0, 0, 0; 0, t49 * t48, 0, 0, 0, 0; t47 * t56 - t52 * t48, -t47 * t53, 0, 0, 0, 0; -t47 * t54 - t50 * t48, -t47 * t55, 0, 0, 0, 0; 0, -t49 * t47, 0, 0, 0, 0; -t55, -t54, 0, 0, 0, 0; t53, -t56, 0, 0, 0, 0; 0, t51, 0, 0, 0, 0;];
JR_rot  = t1;
