% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR7
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
% Datum: 2019-02-26 20:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:52:12
% EndTime: 2019-02-26 20:52:12
% DurationCPUTime: 0.04s
% Computational Cost: add. (80->19), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t95 = qJ(3) + pkin(10) + qJ(5);
t93 = sin(t95);
t98 = cos(qJ(6));
t107 = t93 * t98;
t96 = sin(qJ(6));
t97 = sin(qJ(1));
t106 = t97 * t96;
t105 = t97 * t98;
t99 = cos(qJ(1));
t104 = t99 * t93;
t103 = t99 * t96;
t102 = t99 * t98;
t94 = cos(t95);
t101 = t94 * t106;
t100 = t94 * t102;
t92 = t97 * t93;
t91 = t93 * t96;
t90 = t94 * t103;
t89 = t94 * t105;
t88 = t93 * t102 - t106;
t87 = t93 * t103 + t105;
t86 = t93 * t105 + t103;
t85 = -t93 * t106 + t102;
t1 = [t88, 0, t89, 0, t89, t85; t86, 0, -t100, 0, -t100, t87; 0, 0, -t107, 0, -t107, -t94 * t96; -t87, 0, -t101, 0, -t101, -t86; t85, 0, t90, 0, t90, t88; 0, 0, t91, 0, t91, -t94 * t98; -t99 * t94, 0, t92, 0, t92, 0; -t97 * t94, 0, -t104, 0, -t104, 0; 0, 0, t94, 0, t94, 0;];
JR_rot  = t1;
