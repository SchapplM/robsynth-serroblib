% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:56
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:56:33
% EndTime: 2019-02-22 10:56:33
% DurationCPUTime: 0.04s
% Computational Cost: add. (55->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
t108 = sin(qJ(3));
t109 = sin(qJ(1));
t115 = t109 * t108;
t107 = qJ(4) + qJ(5);
t105 = sin(t107);
t110 = cos(qJ(3));
t114 = t110 * t105;
t106 = cos(t107);
t103 = t110 * t106;
t111 = cos(qJ(1));
t113 = t111 * t108;
t112 = t111 * t110;
t101 = -t109 * t105 + t106 * t113;
t100 = t105 * t113 + t109 * t106;
t99 = t111 * t105 + t106 * t115;
t98 = t105 * t115 - t111 * t106;
t1 = [t101, 0, t109 * t103, -t98, -t98, 0; t99, 0, -t106 * t112, t100, t100, 0; 0, 0, -t108 * t106, -t114, -t114, 0; -t112, 0, t115, 0, 0, 0; -t109 * t110, 0, -t113, 0, 0, 0; 0, 0, t110, 0, 0, 0; t100, 0, t109 * t114, t99, t99, 0; t98, 0, -t105 * t112, -t101, -t101, 0; 0, 0, -t108 * t105, t103, t103, 0;];
JR_rot  = t1;
