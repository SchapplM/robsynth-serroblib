% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:48
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiR_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiR_rot_3_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:47:57
% EndTime: 2019-02-22 11:47:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (38->19), mult. (116->49), div. (0->0), fcn. (167->10), ass. (0->24)
t86 = sin(pkin(14));
t90 = cos(pkin(7));
t104 = t86 * t90;
t88 = sin(pkin(6));
t92 = sin(qJ(1));
t103 = t88 * t92;
t94 = cos(qJ(1));
t102 = t88 * t94;
t89 = cos(pkin(14));
t101 = t89 * t90;
t91 = sin(qJ(2));
t100 = t90 * t91;
t99 = cos(pkin(6));
t98 = t92 * t99;
t97 = t94 * t99;
t93 = cos(qJ(2));
t80 = t92 * t91 - t93 * t97;
t87 = sin(pkin(7));
t96 = t87 * t102 + t80 * t90;
t82 = -t94 * t91 - t93 * t98;
t95 = t87 * t103 + t82 * t90;
t83 = -t91 * t98 + t94 * t93;
t81 = -t91 * t97 - t92 * t93;
t1 = [t81 * t89 + t96 * t86, -t83 * t104 + t82 * t89, 0, 0, 0, 0; t83 * t89 + t95 * t86, t81 * t104 - t80 * t89, 0, 0, 0, 0; 0 (-t86 * t100 + t89 * t93) * t88, 0, 0, 0, 0; -t81 * t86 + t96 * t89, -t83 * t101 - t82 * t86, 0, 0, 0, 0; -t83 * t86 + t95 * t89, t81 * t101 + t80 * t86, 0, 0, 0, 0; 0 (-t89 * t100 - t86 * t93) * t88, 0, 0, 0, 0; t90 * t102 - t80 * t87, t83 * t87, 0, 0, 0, 0; t90 * t103 - t82 * t87, -t81 * t87, 0, 0, 0, 0; 0, t88 * t91 * t87, 0, 0, 0, 0;];
JR_rot  = t1;
