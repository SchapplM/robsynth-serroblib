% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:18
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:18:17
% EndTime: 2019-02-22 12:18:17
% DurationCPUTime: 0.04s
% Computational Cost: add. (80->15), mult. (50->18), div. (0->0), fcn. (87->6), ass. (0->20)
t93 = qJ(2) + qJ(3) + qJ(4);
t92 = cos(t93);
t94 = sin(pkin(11));
t104 = t92 * t94;
t96 = sin(qJ(1));
t103 = t96 * t94;
t95 = cos(pkin(11));
t102 = t96 * t95;
t97 = cos(qJ(1));
t101 = t97 * t94;
t100 = t97 * t95;
t91 = sin(t93);
t99 = t91 * t102;
t98 = t91 * t100;
t90 = t97 * t92;
t89 = t96 * t92;
t88 = t92 * t95;
t87 = t91 * t101;
t86 = t91 * t103;
t1 = [-t92 * t102 + t101, -t98, -t98, -t98, 0, 0; t92 * t100 + t103, -t99, -t99, -t99, 0, 0; 0, t88, t88, t88, 0, 0; t92 * t103 + t100, t87, t87, t87, 0, 0; -t92 * t101 + t102, t86, t86, t86, 0, 0; 0, -t104, -t104, -t104, 0, 0; -t96 * t91, t90, t90, t90, 0, 0; t97 * t91, t89, t89, t89, 0, 0; 0, t91, t91, t91, 0, 0;];
JR_rot  = t1;
