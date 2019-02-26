% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRPRR1_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:48
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.04s
% Computational Cost: add. (24->15), mult. (70->35), div. (0->0), fcn. (100->12), ass. (0->23)
t108 = sin(pkin(11));
t115 = cos(pkin(6));
t122 = t108 * t115;
t109 = sin(pkin(7));
t106 = sin(pkin(13));
t111 = cos(pkin(13));
t116 = sin(qJ(3));
t117 = cos(qJ(3));
t118 = -t106 * t116 + t111 * t117;
t101 = t118 * t109;
t110 = sin(pkin(6));
t121 = t110 * t101;
t114 = cos(pkin(7));
t120 = t110 * t114;
t113 = cos(pkin(11));
t119 = t113 * t115;
t112 = cos(pkin(12));
t107 = sin(pkin(12));
t105 = -t117 * t106 - t116 * t111;
t104 = -t113 * t107 - t112 * t122;
t103 = -t108 * t107 + t112 * t119;
t102 = t118 * t114;
t1 = [0, 0, -t104 * t109 + t108 * t120, 0 -(-t107 * t122 + t113 * t112) * t105 - t104 * t102 - t108 * t121, 0; 0, 0, -t103 * t109 - t113 * t120, 0 -(t107 * t119 + t108 * t112) * t105 - t103 * t102 + t113 * t121, 0; 0, 0, -t110 * t112 * t109 + t115 * t114, 0, -t115 * t101 + (-t102 * t112 - t105 * t107) * t110, 0;];
Jg_rot  = t1;
