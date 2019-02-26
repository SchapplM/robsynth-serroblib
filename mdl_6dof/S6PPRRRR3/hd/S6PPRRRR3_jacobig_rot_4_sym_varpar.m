% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRR3_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobig_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobig_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:05
% EndTime: 2019-02-26 19:44:05
% DurationCPUTime: 0.04s
% Computational Cost: add. (23->16), mult. (67->38), div. (0->0), fcn. (96->12), ass. (0->22)
t102 = sin(pkin(13));
t110 = cos(pkin(6));
t116 = t102 * t110;
t104 = sin(pkin(7));
t105 = sin(pkin(6));
t115 = t104 * t105;
t109 = cos(pkin(7));
t114 = t105 * t109;
t107 = cos(pkin(13));
t113 = t107 * t110;
t112 = cos(qJ(3));
t111 = sin(qJ(3));
t108 = cos(pkin(8));
t106 = cos(pkin(14));
t103 = sin(pkin(8));
t101 = sin(pkin(14));
t100 = -t107 * t101 - t106 * t116;
t99 = -t102 * t101 + t106 * t113;
t98 = -t106 * t115 + t110 * t109;
t97 = -t100 * t104 + t102 * t114;
t96 = -t99 * t104 - t107 * t114;
t1 = [0, 0, t97 -(-(-t101 * t116 + t107 * t106) * t111 + (t100 * t109 + t102 * t115) * t112) * t103 + t97 * t108, 0, 0; 0, 0, t96 -(-(t101 * t113 + t102 * t106) * t111 + (-t107 * t115 + t99 * t109) * t112) * t103 + t96 * t108, 0, 0; 0, 0, t98 -(t110 * t104 * t112 + (t106 * t109 * t112 - t101 * t111) * t105) * t103 + t98 * t108, 0, 0;];
Jg_rot  = t1;
