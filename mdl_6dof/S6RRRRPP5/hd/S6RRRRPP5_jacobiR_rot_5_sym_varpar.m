% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:14
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:14:37
% EndTime: 2019-02-22 12:14:37
% DurationCPUTime: 0.04s
% Computational Cost: add. (53->15), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t105 = qJ(3) + qJ(4);
t103 = sin(t105);
t106 = sin(qJ(2));
t115 = t106 * t103;
t107 = sin(qJ(1));
t114 = t107 * t106;
t108 = cos(qJ(2));
t113 = t108 * t103;
t104 = cos(t105);
t112 = t108 * t104;
t109 = cos(qJ(1));
t111 = t109 * t106;
t110 = t109 * t108;
t101 = t106 * t104;
t100 = t107 * t103 + t104 * t110;
t99 = t103 * t110 - t107 * t104;
t98 = -t109 * t103 + t107 * t112;
t97 = -t109 * t104 - t107 * t113;
t1 = [-t98, -t104 * t111, -t99, -t99, 0, 0; t100, -t104 * t114, t97, t97, 0, 0; 0, t112, -t115, -t115, 0, 0; -t114, t110, 0, 0, 0, 0; t111, t107 * t108, 0, 0, 0, 0; 0, t106, 0, 0, 0, 0; t97, -t103 * t111, t100, t100, 0, 0; t99, -t103 * t114, t98, t98, 0, 0; 0, t113, t101, t101, 0, 0;];
JR_rot  = t1;
