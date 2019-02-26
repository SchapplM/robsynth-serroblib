% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Jg_rot = S6PRRRPP1_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobig_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:44
% EndTime: 2019-02-26 20:08:44
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
t91 = sin(pkin(10));
t92 = sin(pkin(6));
t102 = t91 * t92;
t93 = cos(pkin(10));
t101 = t93 * t92;
t94 = cos(pkin(6));
t96 = sin(qJ(2));
t100 = t94 * t96;
t98 = cos(qJ(2));
t99 = t94 * t98;
t97 = cos(qJ(3));
t95 = sin(qJ(3));
t1 = [0, t102, t91 * t99 + t93 * t96 (-t91 * t100 + t93 * t98) * t95 - t97 * t102, 0, 0; 0, -t101, t91 * t96 - t93 * t99 (t93 * t100 + t91 * t98) * t95 + t97 * t101, 0, 0; 0, t94, -t92 * t98, t92 * t96 * t95 - t94 * t97, 0, 0;];
Jg_rot  = t1;
