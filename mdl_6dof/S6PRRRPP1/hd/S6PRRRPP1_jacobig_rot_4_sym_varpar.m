% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPP1_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobig_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobig_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:54
% EndTime: 2019-02-26 20:08:54
% DurationCPUTime: 0.02s
% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
t83 = sin(pkin(10));
t84 = sin(pkin(6));
t94 = t83 * t84;
t85 = cos(pkin(10));
t93 = t85 * t84;
t86 = cos(pkin(6));
t88 = sin(qJ(2));
t92 = t86 * t88;
t90 = cos(qJ(2));
t91 = t86 * t90;
t89 = cos(qJ(3));
t87 = sin(qJ(3));
t1 = [0, t94, t83 * t91 + t85 * t88 (-t83 * t92 + t85 * t90) * t87 - t89 * t94, 0, 0; 0, -t93, t83 * t88 - t85 * t91 (t83 * t90 + t85 * t92) * t87 + t89 * t93, 0, 0; 0, t86, -t84 * t90, t84 * t88 * t87 - t86 * t89, 0, 0;];
Jg_rot  = t1;
