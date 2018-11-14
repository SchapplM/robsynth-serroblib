% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 0 (0=Basis) von
% S3RRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
%
% Output:
% Jg_rot [3x3]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S3RRR1_jacobig_rot_0_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_jacobig_rot_0_floatb_twist_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_jacobig_rot_0_floatb_twist_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_0_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:16:10
% EndTime: 2018-11-14 10:16:10
% DurationCPUTime: 0.01s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [0, 0, 0; 0, 0, 0; 0, 0, 0;];
Jg_rot  = t1;
