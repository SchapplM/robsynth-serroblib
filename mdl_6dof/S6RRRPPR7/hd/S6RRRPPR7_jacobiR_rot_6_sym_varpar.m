% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:01
% EndTime: 2019-02-26 22:07:01
% DurationCPUTime: 0.12s
% Computational Cost: add. (94->27), mult. (148->32), div. (0->0), fcn. (221->8), ass. (0->24)
t101 = pkin(10) + qJ(6);
t100 = cos(t101);
t102 = sin(qJ(3));
t105 = cos(qJ(3));
t107 = cos(qJ(1));
t104 = sin(qJ(1));
t106 = cos(qJ(2));
t114 = t104 * t106;
t93 = t102 * t114 + t107 * t105;
t94 = -t107 * t102 + t105 * t114;
t99 = sin(t101);
t117 = t93 * t100 - t94 * t99;
t103 = sin(qJ(2));
t110 = t100 * t105 + t102 * t99;
t116 = t110 * t103;
t113 = t107 * t106;
t111 = t94 * t100 + t93 * t99;
t95 = t102 * t113 - t104 * t105;
t96 = t104 * t102 + t105 * t113;
t89 = t96 * t100 + t95 * t99;
t88 = t95 * t100 - t96 * t99;
t109 = t100 * t102 - t105 * t99;
t90 = t109 * t103;
t1 = [-t111, -t107 * t116, -t88, 0, 0, t88; t89, -t104 * t116, -t117, 0, 0, t117; 0, t110 * t106, -t90, 0, 0, t90; -t117, -t107 * t90, t89, 0, 0, -t89; t88, -t104 * t90, t111, 0, 0, -t111; 0, t109 * t106, t116, 0, 0, -t116; t104 * t103, -t113, 0, 0, 0, 0; -t107 * t103, -t114, 0, 0, 0, 0; 0, -t103, 0, 0, 0, 0;];
JR_rot  = t1;
