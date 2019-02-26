% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP6_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobiR_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:32
% EndTime: 2019-02-26 20:46:32
% DurationCPUTime: 0.02s
% Computational Cost: add. (13->4), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
t47 = pkin(9) + qJ(3);
t45 = sin(t47);
t48 = sin(qJ(1));
t51 = t48 * t45;
t46 = cos(t47);
t49 = cos(qJ(1));
t50 = t49 * t46;
t44 = t49 * t45;
t43 = t48 * t46;
t1 = [t49, 0, 0, 0, 0, 0; t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t43, 0, t44, 0, 0, 0; -t50, 0, t51, 0, 0, 0; 0, 0, -t46, 0, 0, 0; -t51, 0, t50, 0, 0, 0; t44, 0, t43, 0, 0, 0; 0, 0, t45, 0, 0, 0;];
JR_rot  = t1;
