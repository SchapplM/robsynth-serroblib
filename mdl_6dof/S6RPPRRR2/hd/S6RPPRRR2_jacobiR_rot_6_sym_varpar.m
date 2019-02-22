% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:18
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:18:19
% EndTime: 2019-02-22 10:18:19
% DurationCPUTime: 0.04s
% Computational Cost: add. (113->19), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->22)
t97 = pkin(11) + qJ(4);
t91 = sin(t97);
t99 = qJ(5) + qJ(6);
t95 = sin(t99);
t107 = t91 * t95;
t96 = cos(t99);
t106 = t91 * t96;
t98 = qJ(1) + pkin(10);
t92 = sin(t98);
t105 = t92 * t95;
t104 = t92 * t96;
t93 = cos(t97);
t103 = t93 * t95;
t102 = t93 * t96;
t94 = cos(t98);
t101 = t94 * t95;
t100 = t94 * t96;
t90 = t93 * t100 + t105;
t89 = -t93 * t101 + t104;
t88 = -t92 * t102 + t101;
t87 = t92 * t103 + t100;
t1 = [t88, 0, 0, -t91 * t100, t89, t89; t90, 0, 0, -t91 * t104, -t87, -t87; 0, 0, 0, t102, -t107, -t107; t87, 0, 0, t91 * t101, -t90, -t90; t89, 0, 0, t91 * t105, t88, t88; 0, 0, 0, -t103, -t106, -t106; -t92 * t91, 0, 0, t94 * t93, 0, 0; t94 * t91, 0, 0, t92 * t93, 0, 0; 0, 0, 0, t91, 0, 0;];
JR_rot  = t1;
