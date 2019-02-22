% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:06
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:06:18
% EndTime: 2019-02-22 11:06:18
% DurationCPUTime: 0.06s
% Computational Cost: add. (73->22), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->26)
t91 = sin(pkin(10));
t92 = cos(pkin(10));
t93 = sin(qJ(6));
t95 = cos(qJ(6));
t100 = t91 * t95 - t92 * t93;
t90 = qJ(2) + pkin(9);
t88 = sin(t90);
t98 = t100 * t88;
t94 = sin(qJ(1));
t106 = t94 * t91;
t105 = t94 * t92;
t96 = cos(qJ(1));
t104 = t96 * t91;
t103 = t96 * t92;
t89 = cos(t90);
t83 = t89 * t106 + t103;
t84 = t89 * t105 - t104;
t102 = t83 * t95 - t84 * t93;
t101 = -t83 * t93 - t84 * t95;
t99 = t91 * t93 + t92 * t95;
t97 = t99 * t88;
t86 = t89 * t103 + t106;
t85 = t89 * t104 - t105;
t82 = t85 * t93 + t86 * t95;
t81 = t85 * t95 - t86 * t93;
t1 = [t101, -t96 * t97, 0, 0, 0, t81; t82, -t94 * t97, 0, 0, 0, t102; 0, t99 * t89, 0, 0, 0, t98; -t102, -t96 * t98, 0, 0, 0, -t82; t81, -t94 * t98, 0, 0, 0, t101; 0, t100 * t89, 0, 0, 0, -t97; t94 * t88, -t96 * t89, 0, 0, 0, 0; -t96 * t88, -t94 * t89, 0, 0, 0, 0; 0, -t88, 0, 0, 0, 0;];
JR_rot  = t1;
