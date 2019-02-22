% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:19
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:19:32
% EndTime: 2019-02-22 11:19:32
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t98 = qJ(2) + pkin(9);
t94 = sin(t98);
t99 = sin(qJ(1));
t106 = t99 * t94;
t97 = qJ(4) + pkin(10);
t95 = cos(t97);
t105 = t99 * t95;
t96 = cos(t98);
t104 = t99 * t96;
t100 = cos(qJ(1));
t103 = t100 * t94;
t102 = t100 * t95;
t101 = t100 * t96;
t93 = sin(t97);
t92 = t95 * t101 + t99 * t93;
t91 = t93 * t101 - t105;
t90 = -t100 * t93 + t95 * t104;
t89 = -t93 * t104 - t102;
t1 = [-t90, -t94 * t102, 0, -t91, 0, 0; t92, -t94 * t105, 0, t89, 0, 0; 0, t96 * t95, 0, -t94 * t93, 0, 0; -t106, t101, 0, 0, 0, 0; t103, t104, 0, 0, 0, 0; 0, t94, 0, 0, 0, 0; t89, -t93 * t103, 0, t92, 0, 0; t91, -t93 * t106, 0, t90, 0, 0; 0, t96 * t93, 0, t94 * t95, 0, 0;];
JR_rot  = t1;
