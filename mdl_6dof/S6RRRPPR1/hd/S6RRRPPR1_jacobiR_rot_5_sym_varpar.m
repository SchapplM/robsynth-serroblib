% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:49
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:49:34
% EndTime: 2019-02-22 11:49:34
% DurationCPUTime: 0.11s
% Computational Cost: add. (59->12), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
t92 = qJ(2) + qJ(3) + pkin(10);
t91 = cos(t92);
t93 = sin(pkin(11));
t103 = t91 * t93;
t95 = sin(qJ(1));
t102 = t95 * t93;
t94 = cos(pkin(11));
t101 = t95 * t94;
t96 = cos(qJ(1));
t100 = t96 * t93;
t99 = t96 * t94;
t90 = sin(t92);
t98 = t90 * t101;
t97 = t90 * t99;
t89 = t96 * t91;
t88 = t95 * t91;
t87 = t91 * t94;
t86 = t90 * t100;
t85 = t90 * t102;
t1 = [-t91 * t101 + t100, -t97, -t97, 0, 0, 0; t91 * t99 + t102, -t98, -t98, 0, 0, 0; 0, t87, t87, 0, 0, 0; t91 * t102 + t99, t86, t86, 0, 0, 0; -t91 * t100 + t101, t85, t85, 0, 0, 0; 0, -t103, -t103, 0, 0, 0; -t95 * t90, t89, t89, 0, 0, 0; t96 * t90, t88, t88, 0, 0, 0; 0, t90, t90, 0, 0, 0;];
JR_rot  = t1;
