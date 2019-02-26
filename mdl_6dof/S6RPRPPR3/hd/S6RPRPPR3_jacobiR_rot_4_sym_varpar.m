% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobiR_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:19
% EndTime: 2019-02-26 20:40:19
% DurationCPUTime: 0.02s
% Computational Cost: add. (14->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
t43 = qJ(1) + pkin(9);
t41 = sin(t43);
t44 = sin(qJ(3));
t48 = t41 * t44;
t45 = cos(qJ(3));
t47 = t41 * t45;
t42 = cos(t43);
t46 = t42 * t44;
t40 = t42 * t45;
t1 = [-t47, 0, -t46, 0, 0, 0; t40, 0, -t48, 0, 0, 0; 0, 0, t45, 0, 0, 0; t42, 0, 0, 0, 0, 0; t41, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t48, 0, t40, 0, 0, 0; t46, 0, t47, 0, 0, 0; 0, 0, t44, 0, 0, 0;];
JR_rot  = t1;
