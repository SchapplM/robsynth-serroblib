% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR14_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:26
% EndTime: 2019-02-26 21:45:26
% DurationCPUTime: 0.04s
% Computational Cost: add. (32->21), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->25)
t100 = sin(pkin(6));
t102 = sin(qJ(4));
t117 = t100 * t102;
t105 = cos(qJ(4));
t116 = t100 * t105;
t106 = cos(qJ(2));
t115 = t100 * t106;
t107 = cos(qJ(1));
t114 = t100 * t107;
t103 = sin(qJ(2));
t104 = sin(qJ(1));
t113 = t104 * t103;
t112 = t104 * t106;
t111 = t107 * t103;
t110 = t107 * t106;
t101 = cos(pkin(6));
t95 = -t101 * t110 + t113;
t109 = t95 * t102 - t105 * t114;
t108 = t102 * t114 + t95 * t105;
t98 = -t101 * t113 + t110;
t97 = t101 * t112 + t111;
t96 = t101 * t111 + t112;
t94 = t97 * t102 + t104 * t116;
t93 = t104 * t117 - t97 * t105;
t1 = [-t96, -t97, 0, 0, 0, 0; t98, -t95, 0, 0, 0, 0; 0, t115, 0, 0, 0, 0; t109, -t98 * t102, 0, t93, 0, 0; -t94, -t96 * t102, 0, -t108, 0, 0; 0, -t103 * t117, 0, t101 * t102 + t105 * t115, 0, 0; t108, -t98 * t105, 0, t94, 0, 0; t93, -t96 * t105, 0, t109, 0, 0; 0, -t103 * t116, 0, t101 * t105 - t102 * t115, 0, 0;];
JR_rot  = t1;
