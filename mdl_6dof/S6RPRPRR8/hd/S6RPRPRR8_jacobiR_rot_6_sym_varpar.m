% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:34
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:34:52
% EndTime: 2019-02-22 10:34:52
% DurationCPUTime: 0.04s
% Computational Cost: add. (83->19), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
t94 = qJ(3) + pkin(10);
t91 = cos(t94);
t95 = qJ(5) + qJ(6);
t92 = sin(t95);
t103 = t91 * t92;
t93 = cos(t95);
t102 = t91 * t93;
t96 = sin(qJ(1));
t101 = t96 * t92;
t100 = t96 * t93;
t97 = cos(qJ(1));
t99 = t97 * t92;
t98 = t97 * t93;
t90 = sin(t94);
t89 = t90 * t98 - t101;
t88 = t90 * t99 + t100;
t87 = t90 * t100 + t99;
t86 = -t90 * t101 + t98;
t1 = [t89, 0, t91 * t100, 0, t86, t86; t87, 0, -t91 * t98, 0, t88, t88; 0, 0, -t90 * t93, 0, -t103, -t103; -t88, 0, -t91 * t101, 0, -t87, -t87; t86, 0, t91 * t99, 0, t89, t89; 0, 0, t90 * t92, 0, -t102, -t102; -t97 * t91, 0, t96 * t90, 0, 0, 0; -t96 * t91, 0, -t97 * t90, 0, 0, 0; 0, 0, t91, 0, 0, 0;];
JR_rot  = t1;
