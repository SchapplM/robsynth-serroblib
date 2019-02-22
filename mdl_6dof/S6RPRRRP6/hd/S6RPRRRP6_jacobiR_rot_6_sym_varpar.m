% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP6
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

function JR_rot = S6RPRRRP6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:53:58
% EndTime: 2019-02-22 10:53:58
% DurationCPUTime: 0.04s
% Computational Cost: add. (81->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t93 = pkin(10) + qJ(3);
t89 = sin(t93);
t94 = qJ(4) + qJ(5);
t91 = sin(t94);
t102 = t89 * t91;
t92 = cos(t94);
t101 = t89 * t92;
t95 = sin(qJ(1));
t100 = t95 * t91;
t99 = t95 * t92;
t96 = cos(qJ(1));
t98 = t96 * t91;
t97 = t96 * t92;
t90 = cos(t93);
t88 = t90 * t97 + t100;
t87 = -t90 * t98 + t99;
t86 = -t90 * t99 + t98;
t85 = t90 * t100 + t97;
t1 = [t86, 0, -t89 * t97, t87, t87, 0; t88, 0, -t89 * t99, -t85, -t85, 0; 0, 0, t90 * t92, -t102, -t102, 0; t85, 0, t89 * t98, -t88, -t88, 0; t87, 0, t89 * t100, t86, t86, 0; 0, 0, -t90 * t91, -t101, -t101, 0; -t95 * t89, 0, t96 * t90, 0, 0, 0; t96 * t89, 0, t95 * t90, 0, 0, 0; 0, 0, t89, 0, 0, 0;];
JR_rot  = t1;
