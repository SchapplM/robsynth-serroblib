% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:55
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobiR_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:55:04
% EndTime: 2019-02-22 09:55:04
% DurationCPUTime: 0.05s
% Computational Cost: add. (59->14), mult. (87->32), div. (0->0), fcn. (134->8), ass. (0->26)
t93 = sin(pkin(11));
t94 = sin(pkin(6));
t104 = t93 * t94;
t95 = cos(pkin(11));
t103 = t94 * t95;
t97 = sin(qJ(2));
t102 = t94 * t97;
t98 = cos(qJ(2));
t101 = t94 * t98;
t96 = cos(pkin(6));
t100 = t96 * t97;
t99 = t96 * t98;
t92 = qJ(3) + qJ(4);
t91 = cos(t92);
t90 = sin(t92);
t89 = -t93 * t100 + t95 * t98;
t88 = -t93 * t99 - t95 * t97;
t87 = t95 * t100 + t93 * t98;
t86 = -t93 * t97 + t95 * t99;
t85 = -t91 * t102 - t96 * t90;
t84 = -t90 * t102 + t96 * t91;
t83 = -t90 * t104 - t89 * t91;
t82 = t91 * t104 - t89 * t90;
t81 = t90 * t103 - t87 * t91;
t80 = -t91 * t103 - t87 * t90;
t1 = [0, t88 * t91, t82, t82, 0, 0; 0, t86 * t91, t80, t80, 0, 0; 0, t91 * t101, t84, t84, 0, 0; 0, -t88 * t90, t83, t83, 0, 0; 0, -t86 * t90, t81, t81, 0, 0; 0, -t90 * t101, t85, t85, 0, 0; 0, t89, 0, 0, 0, 0; 0, t87, 0, 0, 0, 0; 0, t102, 0, 0, 0, 0;];
JR_rot  = t1;
