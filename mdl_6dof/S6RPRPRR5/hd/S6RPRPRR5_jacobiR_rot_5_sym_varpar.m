% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:33
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:32:54
% EndTime: 2019-02-22 10:32:54
% DurationCPUTime: 0.04s
% Computational Cost: add. (50->13), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->14)
t41 = pkin(10) + qJ(3);
t39 = sin(t41);
t40 = cos(t41);
t42 = sin(qJ(5));
t44 = cos(qJ(5));
t49 = -t39 * t44 + t40 * t42;
t46 = t39 * t42 + t40 * t44;
t45 = cos(qJ(1));
t43 = sin(qJ(1));
t35 = t46 * t45;
t34 = t49 * t45;
t33 = t46 * t43;
t32 = t49 * t43;
t1 = [-t33, 0, t34, 0, -t34, 0; t35, 0, t32, 0, -t32, 0; 0, 0, t46, 0, -t46, 0; t32, 0, t35, 0, -t35, 0; -t34, 0, t33, 0, -t33, 0; 0, 0, -t49, 0, t49, 0; -t45, 0, 0, 0, 0, 0; -t43, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
