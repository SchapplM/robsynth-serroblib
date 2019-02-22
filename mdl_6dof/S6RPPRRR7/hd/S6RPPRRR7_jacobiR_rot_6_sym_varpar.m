% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:21
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:20:57
% EndTime: 2019-02-22 10:20:57
% DurationCPUTime: 0.04s
% Computational Cost: add. (80->19), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t96 = pkin(10) + qJ(4) + qJ(5);
t94 = sin(t96);
t99 = cos(qJ(6));
t108 = t94 * t99;
t97 = sin(qJ(6));
t98 = sin(qJ(1));
t107 = t98 * t97;
t106 = t98 * t99;
t100 = cos(qJ(1));
t105 = t100 * t94;
t104 = t100 * t97;
t103 = t100 * t99;
t95 = cos(t96);
t102 = t95 * t107;
t101 = t95 * t103;
t93 = t98 * t94;
t92 = t94 * t97;
t91 = t95 * t104;
t90 = t95 * t106;
t89 = t94 * t103 - t107;
t88 = t94 * t104 + t106;
t87 = t94 * t106 + t104;
t86 = -t94 * t107 + t103;
t1 = [t89, 0, 0, t90, t90, t86; t87, 0, 0, -t101, -t101, t88; 0, 0, 0, -t108, -t108, -t95 * t97; -t88, 0, 0, -t102, -t102, -t87; t86, 0, 0, t91, t91, t89; 0, 0, 0, t92, t92, -t95 * t99; -t100 * t95, 0, 0, t93, t93, 0; -t98 * t95, 0, 0, -t105, -t105, 0; 0, 0, 0, t95, t95, 0;];
JR_rot  = t1;
