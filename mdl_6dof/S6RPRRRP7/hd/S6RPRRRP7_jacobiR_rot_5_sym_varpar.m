% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:54
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:54:35
% EndTime: 2019-02-22 10:54:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (81->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t92 = pkin(10) + qJ(3);
t88 = sin(t92);
t93 = qJ(4) + qJ(5);
t90 = sin(t93);
t101 = t88 * t90;
t91 = cos(t93);
t100 = t88 * t91;
t94 = sin(qJ(1));
t99 = t94 * t90;
t98 = t94 * t91;
t95 = cos(qJ(1));
t97 = t95 * t90;
t96 = t95 * t91;
t89 = cos(t92);
t87 = t89 * t96 + t99;
t86 = -t89 * t97 + t98;
t85 = -t89 * t98 + t97;
t84 = t89 * t99 + t96;
t1 = [t85, 0, -t88 * t96, t86, t86, 0; t87, 0, -t88 * t98, -t84, -t84, 0; 0, 0, t89 * t91, -t101, -t101, 0; t84, 0, t88 * t97, -t87, -t87, 0; t86, 0, t88 * t99, t85, t85, 0; 0, 0, -t89 * t90, -t100, -t100, 0; -t94 * t88, 0, t95 * t89, 0, 0, 0; t95 * t88, 0, t94 * t89, 0, 0, 0; 0, 0, t88, 0, 0, 0;];
JR_rot  = t1;
