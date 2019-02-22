% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:17
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR9_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:17:52
% EndTime: 2019-02-22 11:17:52
% DurationCPUTime: 0.05s
% Computational Cost: add. (27->17), mult. (87->31), div. (0->0), fcn. (134->8), ass. (0->25)
t86 = sin(pkin(6));
t89 = sin(qJ(2));
t103 = t86 * t89;
t91 = cos(qJ(5));
t102 = t86 * t91;
t92 = cos(qJ(2));
t101 = t86 * t92;
t93 = cos(qJ(1));
t100 = t86 * t93;
t90 = sin(qJ(1));
t99 = t90 * t89;
t98 = t90 * t92;
t97 = t93 * t89;
t96 = t93 * t92;
t87 = cos(pkin(6));
t81 = t87 * t97 + t98;
t88 = sin(qJ(5));
t95 = t91 * t100 - t81 * t88;
t94 = t88 * t100 + t81 * t91;
t83 = -t87 * t99 + t96;
t82 = -t87 * t98 - t97;
t80 = -t87 * t96 + t99;
t79 = t90 * t102 + t83 * t88;
t78 = -t90 * t86 * t88 + t83 * t91;
t1 = [t95, t82 * t88, 0, 0, t78, 0; t79, -t80 * t88, 0, 0, t94, 0; 0, t88 * t101, 0, 0, t89 * t102 - t87 * t88, 0; -t94, t82 * t91, 0, 0, -t79, 0; t78, -t80 * t91, 0, 0, t95, 0; 0, t91 * t101, 0, 0, -t88 * t103 - t87 * t91, 0; t80, -t83, 0, 0, 0, 0; t82, -t81, 0, 0, 0, 0; 0, -t103, 0, 0, 0, 0;];
JR_rot  = t1;
