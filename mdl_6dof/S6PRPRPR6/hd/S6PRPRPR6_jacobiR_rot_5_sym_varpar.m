% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:12
% EndTime: 2019-02-26 19:49:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (37->22), mult. (107->51), div. (0->0), fcn. (156->10), ass. (0->26)
t102 = sin(qJ(4));
t96 = sin(pkin(11));
t113 = t102 * t96;
t98 = sin(pkin(6));
t112 = t102 * t98;
t99 = cos(pkin(11));
t111 = t102 * t99;
t104 = cos(qJ(4));
t110 = t104 * t98;
t105 = cos(qJ(2));
t109 = t105 * t98;
t101 = cos(pkin(6));
t103 = sin(qJ(2));
t108 = t101 * t103;
t107 = t101 * t105;
t106 = t102 * t103;
t100 = cos(pkin(10));
t97 = sin(pkin(10));
t94 = -t101 * t102 - t104 * t109;
t93 = t100 * t105 - t97 * t108;
t92 = t100 * t103 + t97 * t107;
t91 = t100 * t108 + t97 * t105;
t90 = -t100 * t107 + t97 * t103;
t89 = t100 * t112 + t90 * t104;
t88 = t92 * t104 - t97 * t112;
t1 = [0, t93 * t111 - t92 * t96, 0, t88 * t99, 0, 0; 0, t91 * t111 - t90 * t96, 0, t89 * t99, 0, 0; 0 (t105 * t96 + t99 * t106) * t98, 0, t94 * t99, 0, 0; 0, -t93 * t113 - t92 * t99, 0, -t88 * t96, 0, 0; 0, -t91 * t113 - t90 * t99, 0, -t89 * t96, 0, 0; 0 (t105 * t99 - t96 * t106) * t98, 0, -t94 * t96, 0, 0; 0, -t93 * t104, 0, t92 * t102 + t97 * t110, 0, 0; 0, -t91 * t104, 0, -t100 * t110 + t90 * t102, 0, 0; 0, -t103 * t110, 0, t101 * t104 - t102 * t109, 0, 0;];
JR_rot  = t1;
