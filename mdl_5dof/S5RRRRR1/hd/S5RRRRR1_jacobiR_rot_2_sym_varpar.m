% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR1_jacobiR_rot_2_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_jacobiR_rot_2_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_jacobiR_rot_2_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:37:45
% EndTime: 2019-02-26 19:37:45
% DurationCPUTime: 0.02s
% Computational Cost: add. (9->9), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
t8 = sin(qJ(2));
t9 = sin(qJ(1));
t15 = t9 * t8;
t11 = cos(qJ(1));
t14 = t11 * t8;
t10 = cos(qJ(2));
t13 = t9 * t10;
t12 = t11 * t10;
t1 = [-t13, -t14, 0, 0, 0; t12, -t15, 0, 0, 0; 0, -t10, 0, 0, 0; t15, -t12, 0, 0, 0; -t14, -t13, 0, 0, 0; 0, t8, 0, 0, 0; -t11, 0, 0, 0, 0; -t9, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
JR_rot  = t1;
