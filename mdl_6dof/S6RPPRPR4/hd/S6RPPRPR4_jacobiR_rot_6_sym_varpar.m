% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:10
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRPR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:10:27
% EndTime: 2019-02-22 10:10:27
% DurationCPUTime: 0.04s
% Computational Cost: add. (57->16), mult. (88->24), div. (0->0), fcn. (141->8), ass. (0->20)
t105 = cos(qJ(1));
t104 = sin(qJ(1));
t93 = qJ(4) + pkin(10);
t91 = sin(t93);
t94 = sin(qJ(6));
t103 = t91 * t94;
t95 = cos(qJ(6));
t102 = t91 * t95;
t92 = cos(t93);
t101 = t92 * t94;
t100 = t92 * t95;
t99 = cos(pkin(9));
t98 = sin(pkin(9));
t86 = -t104 * t98 - t105 * t99;
t87 = -t104 * t99 + t105 * t98;
t97 = t87 * t100 + t86 * t94;
t96 = t87 * t101 - t86 * t95;
t85 = -t86 * t100 + t87 * t94;
t84 = t86 * t101 + t87 * t95;
t1 = [t97, 0, 0, t86 * t102, 0, t84; t85, 0, 0, t87 * t102, 0, t96; 0, 0, 0, -t100, 0, t103; -t96, 0, 0, -t86 * t103, 0, -t85; t84, 0, 0, -t87 * t103, 0, t97; 0, 0, 0, t101, 0, t102; t87 * t91, 0, 0, -t86 * t92, 0, 0; -t86 * t91, 0, 0, -t87 * t92, 0, 0; 0, 0, 0, -t91, 0, 0;];
JR_rot  = t1;
