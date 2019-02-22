% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:09
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:09:14
% EndTime: 2019-02-22 11:09:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t95 = qJ(2) + pkin(9);
t91 = sin(t95);
t96 = sin(qJ(1));
t103 = t96 * t91;
t94 = pkin(10) + qJ(5);
t92 = cos(t94);
t102 = t96 * t92;
t93 = cos(t95);
t101 = t96 * t93;
t97 = cos(qJ(1));
t100 = t97 * t91;
t99 = t97 * t92;
t98 = t97 * t93;
t90 = sin(t94);
t89 = t96 * t90 + t92 * t98;
t88 = t90 * t98 - t102;
t87 = t92 * t101 - t97 * t90;
t86 = -t90 * t101 - t99;
t1 = [-t87, -t91 * t99, 0, 0, -t88, 0; t89, -t91 * t102, 0, 0, t86, 0; 0, t93 * t92, 0, 0, -t91 * t90, 0; -t103, t98, 0, 0, 0, 0; t100, t101, 0, 0, 0, 0; 0, t91, 0, 0, 0, 0; t86, -t90 * t100, 0, 0, t89, 0; t88, -t90 * t103, 0, 0, t87, 0; 0, t93 * t90, 0, 0, t91 * t92, 0;];
JR_rot  = t1;
